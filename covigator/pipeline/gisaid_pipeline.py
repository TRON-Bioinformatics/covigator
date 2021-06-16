#!/usr/bin/env python
import base64
import os
import pathlib
import zlib
from dataclasses import dataclass

from Bio import pairwise2
from Bio import Align
from Bio.pairwise2 import format_alignment
from logzero import logger

from covigator.configuration import Configuration
from covigator.database.model import SampleGisaid, Gene
from covigator.database.queries import Queries

CHROMOSOME = "MN908947.3"


@dataclass
class Variant:
    gene_name: str
    position: int
    reference: str
    alternate: str

    def to_vcf_line(self):
        return self.gene_name, str(self.position), ".", self.reference, self.alternate, "255", ".", ".", "."


class GisaidPipeline:

    MAPPING_GENES = {
        "Spike": "S",
        "M": "M",
        "N": "N",
        "E": "E",
        "MN908947.3": "MN908947.3"
    }

    def __init__(self, config: Configuration, queries: Queries):
        self.config = config
        self.genes = {g.name: g for g in queries.get_genes()}

    def run(self, sample: SampleGisaid):
        logger.info("Processing {}".format(sample.run_accession))
        mutations = []
        for g, s in sample.sequence.items():
            gene = self.genes.get(g)
            # the sequences corresponding to protein domains are not called right now
            if gene is not None:
                alignment = self._run_alignment(
                    sequence=self.decompress_sequence(s).strip("*"), reference=gene.sequence)
                mutations.extend(self._call_mutations(gene=gene, alignment=alignment))
        local_folder = os.path.join(self.config.storage_folder, sample.run_accession)
        if not os.path.exists(local_folder):
            pathlib.Path(local_folder).mkdir(parents=True, exist_ok=True)
        output_vcf = os.path.join(local_folder, "gisaid.vcf")
        
        self._output_vcf(mutations, output_vcf)
        return output_vcf

    @staticmethod
    def decompress_sequence(s):
        return zlib.decompress(base64.b64decode(s)).decode('utf-8')

    def _run_alignment(self, sequence: str, reference: str):
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.match = 2
        aligner.mismatch = -1
        aligner.open_gap_score = -3
        aligner.extend_gap_score = -0.1
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0
        alignments = aligner.align(reference, sequence)
        #alignments = None
        #try:
        #    alignments = pairwise2.align.globalms(sequence, reference, 2, -1, -3, -0.1, one_alignment_only=True)
        #except:
        #    print("Not working")
        aln = alignments[0]
        return aln

    def _call_mutations(self, gene: Gene, alignment):
        #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  /home/sorn/development/ngs_pipelines/covigator/covigator/processor/test_folder/Aligned.out.bam
        #MN908947.3      9924    .       C       T       228     .       DP=139;VDB=0.784386;SGB=-0.693147;RPB=0.696296;MQB=1;MQSB=1;BQB=0.740741;MQ0F=0;AC=1;AN=1;DP4=2,0,123,12;MQ=60  GT:PL   1:255,0
        mutations = []
        alternate = alignment.query
        reference = alignment.target

        opened_insertion = None
        opened_deletion = None
        for i in range(len(reference)):
            if alternate[i] != reference[i]:
                if reference[i] == "-":
                    # insertion
                    if opened_insertion is None:
                        # stores for later
                        opened_insertion = Variant(
                            gene_name=gene.name, position=i, reference="-", alternate=alternate[i])
                    else:
                        # extends
                        opened_insertion.alternate = opened_insertion.alternate + alternate[i]
                elif alternate[i] == "-":
                    # deletion
                    if opened_deletion is None:
                        # stores for later
                        opened_deletion = Variant(
                            gene_name=gene.name, position=i, reference=reference[i], alternate="-")
                    else:
                        # extends
                        opened_deletion.reference = opened_deletion.reference + reference[i]
                else:
                    # substitution
                    mutations.append(Variant(
                        gene_name=gene.name, position=i, reference=reference[i], alternate=alternate[i]))
            else:
                if opened_insertion is not None:
                    mutations.append(opened_insertion)
                    opened_insertion = None
                if opened_deletion is not None:
                    mutations.append(opened_deletion)
                    opened_deletion = None

        # stores also deletion at the very end
        if opened_insertion is not None:
            mutations.append(opened_insertion)
        if opened_deletion is not None:
            mutations.append(opened_deletion)

        return mutations
    
    def _output_vcf(self, mutations, output_vcf):
        with open(output_vcf, "w") as vcf_out:
            header = (
                "##fileformat=VCFv4.0",
                "##FILTER=<ID=PASS,Description=\"All filters passed\">",
                "##contig=<ID=MN908947.3>",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
            )
            for row in header:
                vcf_out.write(row + "\n")
            for row in mutations:
                vcf_out.write("\t".join(row.to_vcf_line()) + "\n")
