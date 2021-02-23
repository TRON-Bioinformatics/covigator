#!/bin/bash
from Bio import SeqIO
from Bio.SeqUtils import IUPACData
from cyvcf2 import VCF, Writer, Variant


def create_variant(chromosome, position, reference_base, alternate_base):
    variant = Variant()
    variant.CHROM = chromosome
    variant.POS = position
    variant.REF = reference_base
    variant.ALT = alternate_base
    return variant


def main():
    output = Writer("Sars_cov_2.all_variants.vcf", VCF("vcf_template.vcf"))
    for record in SeqIO.parse("Sars_cov_2.ASM985889v3.dna.toplevel.fa", "fasta"):
        chromosome = record.name
        position = 0
        for reference_base in record.seq:
            position += 1
            if reference_base not in IUPACData.unambiguous_dna_letters:
                continue    # we don't want mutations on ambiguous bases
            for alternate_base in IUPACData.unambiguous_dna_letters:
                if alternate_base != reference_base:
                    # creates one variant for each alternative variant
                    #output.write_record(
                    #    create_variant(
                    #        chromosome=chromosome, position=position, reference_base=reference_base,
                    #        alternate_base=alternate_base))
                    variant = output.variant_from_string("{}\t{}\t.\t{}\t{}\t.\tPASS\t.".format(
                        chromosome, position, reference_base, alternate_base))
                    output.write_record(variant)

    output.close()


if __name__ == '__main__':
    main()
