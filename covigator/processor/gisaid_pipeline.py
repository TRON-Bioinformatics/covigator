#!/usr/bin/env python

import operator
import sys

from Bio import SeqIO, pairwise2
from Bio.Align import AlignInfo

# TODO: Implement pair-wise alignment with BioPython

class Pipeline:
    
    reference_file = "/projects/SARS-CoV-2/RefSeq_surface_glycoprotein.fasta"
    gisaid_file = "/scratch/info/projects/SARS-CoV-2/gisaid/gisaid_cov2020_sequences_01APR2020_human_host_low_cov_excl.fasta"

    def get_sequence(self, sequence_id):
        for record in SeqIO.parse(self.gisaid_file, "fasta"):
            if sequence_id == record.id:
                return record.seq
        return None

    def run_alignment(self, sequence_id):
        sequence = self.get_sequence(sequence_id)
        reference = SeqIO.read(self.reference_file, "fasta")
        alignments = pairwise2.align.globalxx(sequence, reference.seq)
        #for aln in alignments:
        #    print(aln)
        aln = alignments[0]
        seq1 = aln[0]
        seq2 = aln[1]
        print(aln)

        #print(pairwise2.format_alignment(*alignments[0]))
        return aln


    def call_mutations(self, alignment):
        #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  /home/sorn/development/ngs_pipelines/covigator/covigator/processor/test_folder/Aligned.out.bam
        #MN908947.3      9924    .       C       T       228     .       DP=139;VDB=0.784386;SGB=-0.693147;RPB=0.696296;MQB=1;MQSB=1;BQB=0.740741;MQ0F=0;AC=1;AN=1;DP4=2,0,123,12;MQ=60  GT:PL   1:255,0
        mutations = []
        seq1 = alignment[0]
        seq2 = alignment[1]
        
        for i in range(len(seq2)):
            #print(seq1[i], seq2[i])
            if seq1[i] != seq2[i]:
                mutations.append(("MN908947.3", str(i), ".", seq2[i], seq1[i], "255", ".", ".", "."))
            #vcf_out.write("MN908947.3\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(i, ".", seq2[i], seq1[i], 0, ".", ".", "."))
            #print("MN908947.3\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(i, ".", seq2[i], seq1[i], 0, ".", ".", "."))
                #print(i)
                #print("Alt:", seq1[i])
                #print("Ref:", seq2[i])
        #aln_info = AlignInfo.SummaryInfo(alignment)
        #print(aln_info.dumb_consensus())
        return mutations
    
    def output_vcf(self, mutations):
        with open("gisaid.vcf", "w") as vcf_out:
            vcf_out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")
            for row in mutations:
                vcf_out.write("\t".join(row) + "\n")

    def run(self, run_accession: str):
        print("Processing {}".format(run_accession))
        alignment = run_alignment(run_accession)
        mutations = call_mutations(alignment)
        vcf_file = output_vcf(mutations)

        return vcf_file
                

    #Processor.load()
    #VcfLoader().load(vcf_file=job.vcf_path, sample=job.run_accession, session=session)
    #return vcf_file

if __name__ == "__main__":
    # For now just read in the whole GISAID protein sequences
    # in order to test this
    # Populate DB later and only run for new sequences


    # TODO: Use nucleotide sequences
    #reference_file = "/projects/SARS-CoV-2/RefSeq_surface_glycoprotein_AA.fasta"
    #reference_file = "/projects/SARS-CoV-2/RefSeq_surface_glycoprotein.fasta"
    #reference = SeqIO.read(reference_file, "fasta")
    #print(reference)
    gisaid_file = "/scratch/info/projects/SARS-CoV-2/gisaid/gisaid_cov2020_sequences_01APR2020_human_host_low_cov_excl.fasta"
    sequence_id = "hCoV-19/USA/WA-S95/2020|EPI_ISL_417148|2020-02-28"
    pipe = Pipeline()

    vcf_file = pipe.run(sequence_id)
    #for record in SeqIO.parse(gisaid_file, "fasta"):
    #for record in SeqIO.parse("/projects/SARS-CoV-2/gisaid_cov2020_sequences_01APR2020_human_host_low_cov_excl.S_prot_AA.fasta", "fasta"):

        #vcf_file = pipe.run(record.id, record.seq, reference.seq)
        #break
