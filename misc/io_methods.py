#!/usr/bin/env python

import os
import pwd
import grp
import stat
import pandas as pd

def convert_bed_line_to_exons(words):
    """Returns a list of exons (chr,start,stop,transcript name,exon number,strand). Input is a list representing an eleven field bed line. The start and stop coordinates are returned as zero base closed interval."""
    exon_list = []
    start = int(words[1])
#    cdef unsigned long int end = int(words[2])
    name = words[3]
#    (start,end,name) = words[1:4]
#    start = int(start)
    strand = words[5]
#    print(strand, words[5])
    num_exons = int(words[9])
    exon_starts = words[11].strip(",").split(",")
    exon_sizes = words[10].strip(",").split(",")
    exon_startsI = [start + int(x) for x in exon_starts]
    exon_sizesI = [int(x) for x in exon_sizes]
    r = range(num_exons)
    if strand == "-":
        r = range(num_exons-1,-1,-1)
    for ne in r:
        abs_e_end = exon_startsI[ne] + exon_sizesI[ne] - 1
        exon_number_correct_order = (num_exons - 1 - ne) if strand == "-" else ne
        exon_words = [exon_startsI[ne],abs_e_end,name,exon_number_correct_order,strand]
        exon_list.append(exon_words)
    return exon_list

class IOMethods(object):
    @staticmethod
    def create_folder(path):
        '''This function creates a specified folder, if not already existing and grants the respective permission.'''
        if not os.path.exists(path):
            print("Creating folder", path)
            os.makedirs(path)

    @staticmethod
    def get_fastq_files(paths):
        '''This function returns a list of fastq files for the given list of paths.'''
        fastqs = []
        for path in paths:
            if os.path.isdir(path):
                files = os.listdir(path)
                for file in files:
                    file_path = os.path.join(path,file)
                    if os.path.isfile(file_path) and file.endswith((".fastq.gz","fastq")):
                        fastqs.append(file_path)
                    elif os.path.isdir(file_path):
                        fastqs_tmp = IOMethods.get_fastq_files([file_path])
                        fastqs.extend(fastqs_tmp)
            elif os.path.isfile(path) and path.endswith((".fastq.gz","fastq")):
                fastqs.append(path)
        return fastqs


    @staticmethod
    def load_reference_transcripts(bedfile, gene_symbols_file = None):
        gene_symbols = None
        if gene_symbols_file:
            gene_symbols = IOMethods.load_gene_symbols(gene_symbols_file)

        exons = {}
        transcripts = {}
#        exlength = {}
#        translength = {}
        numt = 0
        nume = 0

        # read bed file
        with open(bedfile, "r") as f:
            for line in f:
                elements = line.rstrip().split("\t")
                chrom = elements[0]
                t_start = int(elements[1])
                t_end = int(elements[2])
                t_strand = elements[5]
                t_id = elements[3]

                if gene_symbols_file:
                    t_id = gene_symbols[t_id]
#                t_ele = (t_start, t_end, t_id)
#                t = Transcript(t_start, t_end, t_id)
#                if chrom in transcripts:
                try:
                    transcripts[chrom].append((t_start, t_end, t_id))
#                    transcripts[chrom].append(t_ele)
#                    transcripts[chrom].append(t)
                except:
#                else:
                    transcripts[chrom] = [(t_start, t_end, t_id)]
#                    transcripts[chrom] = [t_ele]
#                    transcripts[chrom] = [t]
                # convert to exons
                exon_list = convert_bed_line_to_exons(elements)
                for exon in exon_list:
#                    e_start, e_end, e_id, e_num, e_strand = exon
                    e_start = exon[0]
                    e_end = exon[1]
#                    e_id = exon[2]
#                    e_num = exon[3]
#                    e_strand = exon[4]
#                    e = Transcript(e_start, e_end, t_id)
                    e_ele = (e_start, e_end, t_strand, t_id)
                    try:
#                        if e in exons[chrom]
#                        exons[chrom].append((e_start, e_end, t_id))
                        exons[chrom].append(e_ele)
#                        exons[chrom].append(e)
                    except:
#                        exons[chrom] = [(e_start, e_end, t_id)]
                        exons[chrom] = [e_ele]
#                        exons[chrom] = [e]
#        nume = numt = 0
#        for i in exons.keys():

        for i in exons:
#            print(exons[i])
            exons[i] = list(set(exons[i]))
            
#            exons[i] = set(exons[i])
#            exons[i].sort(key = lambda x: x.start, reverse=False)
            exons[i].sort()
            nume += len(exons[i])
#        for i in transcripts.keys():
        for i in transcripts:
#            transcripts[i].sort(key = lamdba x: x.start, reverse=False)
            transcripts[i].sort()
            numt += len(transcripts[i])

        return transcripts, exons, numt, nume

    @staticmethod
    def load_gene_symbols(infile):
        symbol_dict = {}
        with open(infile) as inf:
            for line in inf:
                elements = line.rstrip().split("\t")
                trans_id = elements[0]
                gene_symbol = elements[1]
                        
                symbol_dict[trans_id] = gene_symbol
        return symbol_dict

    @staticmethod
    def load_gene_symbols_pandas(infile):
        pandas_tab = pd.read_table(infile, header=None, names = ["ID", "Gene Symbol", "Length"], usecols=["ID", "Gene Symbol"], sep="\t", index_col=0)
        symbol_dict = pandas_tab.to_dict()
#        print(symbol_dict)
        return symbol_dict


def main():
    pass

if __name__ == '__main__':
    main()
