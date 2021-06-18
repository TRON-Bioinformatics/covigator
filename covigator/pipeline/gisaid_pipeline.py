#!/usr/bin/env python
import os
import pathlib
from dataclasses import dataclass
from Bio import Align
from Bio.Align import PairwiseAlignment
from logzero import logger
from typing import List
from covigator.configuration import Configuration
from covigator.database.model import SampleGisaid
from covigator.database.queries import Queries
from covigator.misc.compression import decompress_sequence

CHROMOSOME = "MN908947.3"


@dataclass
class Variant:
    position: int
    reference: str
    alternate: str

    def to_vcf_line(self):
        # transform 0-based position to 1-based position
        return CHROMOSOME, str(self.position + 1), ".", self.reference, self.alternate, ".", "PASS", "."


class AssemblyVariantCaller:

    def call_variants(self, sequence: str, reference: str) -> List[Variant]:
        alignment = self._run_alignment(sequence=sequence, reference=reference)
        variants = self._call_mutations(alignment)
        return variants

    def _run_alignment(self, sequence: str, reference: str) -> PairwiseAlignment:
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.match = 2
        aligner.mismatch = -1
        aligner.open_gap_score = -3
        aligner.extend_gap_score = -0.1
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0
        alignments = aligner.align(reference, sequence)
        return alignments[0]

    def _call_mutations(self, alignment: PairwiseAlignment) -> List[Variant]:
        # CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
        # MN908947.3      9924    .       C       T       228     .
        # DP=139;VDB=0.784386;SGB=-0.693147;RPB=0.696296;MQB=1;MQSB=1;BQB=0.740741;MQ0F=0;AC=1;AN=1;DP4=2,0,123,12;MQ=60
        # GT:PL   1:255,0
        alternate = alignment.query
        reference = alignment.target

        variants = []
        prev_ref_end = None
        prev_alt_end = None
        for (ref_start, ref_end), (alt_start, alt_end) in zip(alignment.aligned[0], alignment.aligned[1]):
            # calls indels
            # NOTE: it does not call indels at beginning and end of sequence
            if prev_ref_end is not None and prev_ref_end != ref_start:
                # deletion
                if ref_start - prev_ref_end <= 50:  # skips deletions longer than 50 bp
                    ref = reference[prev_ref_end:ref_start + 1]
                    if 'N' not in ref:  # do not call deletions with Ns
                        variants.append(Variant(
                            position=prev_ref_end,
                            reference=ref,
                            alternate=reference[prev_ref_end]))
            elif prev_ref_end is not None and  prev_alt_end != alt_start:
                # insertion
                if alt_start - prev_alt_end <= 50:  # skips insertions longer than 50 bp
                    ref = reference[prev_ref_end]
                    alt = alternate[prev_alt_end:alt_start]
                    if ref != 'N' and 'N' not in alt:   # do not call insertions with Ns
                        variants.append(Variant(
                            position=prev_ref_end,
                            reference=ref,
                            alternate=ref + alt))

            # calls SNVs
            for pos, ref, alt in zip(
                    range(ref_start, ref_end), reference[ref_start: ref_end], alternate[alt_start: alt_end]):
                # contiguous SNVs are reported separately
                if ref != alt and ref != 'N' and alt != 'N':    # do not call SNVs on Ns
                    variants.append(Variant(position=pos, reference=ref, alternate=alt))

            prev_ref_end = ref_end
            prev_alt_end = alt_end

        return variants


class GisaidPipeline:

    def __init__(self, config: Configuration, queries: Queries):
        self.config = config
        self.genes = {g.name: g for g in queries.get_genes()}
        self.variant_caller = AssemblyVariantCaller()

    def run(self, sample: SampleGisaid):
        logger.info("Processing {}".format(sample.run_accession))
        mutations = []
        for g, s in sample.sequence.items():
            gene = self.genes.get(g)
            # the sequences corresponding to protein domains are not called right now
            if gene is not None:
                mutations.extend(self.variant_caller.call_variants(
                    sequence=decompress_sequence(s).strip("*"), reference=gene.sequence))
        local_folder = os.path.join(self.config.storage_folder, sample.run_accession)
        if not os.path.exists(local_folder):
            pathlib.Path(local_folder).mkdir(parents=True, exist_ok=True)
        output_vcf = os.path.join(local_folder, "gisaid.vcf")
        
        self._output_vcf(mutations, output_vcf)
        return output_vcf

    def _output_vcf(self, mutations, output_vcf):
        with open(output_vcf, "w") as vcf_out:
            header = (
                "##fileformat=VCFv4.0",
                "##FILTER=<ID=PASS,Description=\"All filters passed\">",
                "##contig=<ID=MN908947.3>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
            )
            for row in header:
                vcf_out.write(row + "\n")
            for row in mutations:
                vcf_out.write("\t".join(row.to_vcf_line()) + "\n")
