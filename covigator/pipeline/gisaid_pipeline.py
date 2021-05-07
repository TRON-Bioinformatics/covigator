#!/usr/bin/env python
import base64
import os
import zlib

from Bio import SeqIO, pairwise2
from covigator.configuration import Configuration
from covigator.database.model import SampleGisaid


class GisaidPipeline:

    def __init__(self, config: Configuration):
        self.config = config

    def run(self, sample: SampleGisaid):
        print("Processing {}".format(sample.run_accession))
        mutations = []
        for s in sample.sequence.values():
            alignment = self._run_alignment(self.decompress_sequence(s))
            mutations.extend(self._call_mutations(alignment))
        local_folder = os.path.join(self.config.storage_folder, sample.run_accession)
        output_vcf = os.path.join(local_folder, "gisaid.vcf")
        self._output_vcf(mutations, output_vcf)
        return output_vcf

    @staticmethod
    def decompress_sequence(s):
        return zlib.decompress(base64.b64decode(s)).decode('utf-8')

    def _run_alignment(self, sequence):
        reference = SeqIO.read(self.config.reference_genome, "fasta")
        alignments = pairwise2.align.globalxx(sequence, reference.seq)
        aln = alignments[0]
        return aln

    def _call_mutations(self, alignment):
        #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  /home/sorn/development/ngs_pipelines/covigator/covigator/processor/test_folder/Aligned.out.bam
        #MN908947.3      9924    .       C       T       228     .       DP=139;VDB=0.784386;SGB=-0.693147;RPB=0.696296;MQB=1;MQSB=1;BQB=0.740741;MQ0F=0;AC=1;AN=1;DP4=2,0,123,12;MQ=60  GT:PL   1:255,0
        mutations = []
        seq1 = alignment[0]
        seq2 = alignment[1]
        
        for i in range(len(seq2)):
            if seq1[i] != seq2[i]:
                mutations.append(("MN908947.3", str(i), ".", seq2[i], seq1[i], "255", ".", ".", "."))

        return mutations
    
    def _output_vcf(self, mutations, output_vcf):
        with open(output_vcf, "w") as vcf_out:
            vcf_out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")
            for row in mutations:
                vcf_out.write("\t".join(row) + "\n")
