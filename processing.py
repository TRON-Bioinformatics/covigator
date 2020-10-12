#!/usr/bin/env python

# Scripts
# Minion scripts for corona virus pipeline
# /projects/SARS-CoV-2/minion_gridion.sh

# RNA-Seq script for corona virus pipeline
# /projects/SARS-CoV-2/documentation_THBU.sh


# Data
# RNA-Seq data for corona virus pipeline
# /scratch/info/projects/SARS-CoV-2/Novoalign/fastqs/

# Minion data for corona virus pipeline
# /scratch/info/projects/SARS-CoV-2/sra2

# GISAID assembly data for Coronavirus
# /scratch/info/projects/SARS-CoV-2/gisaid/


import config as cfg

from misc.io_methods import IOMethods
from misc.queue import Queue

def download_virus_db(db_path):
    cmds = []
    print("Downloading fasta files.")
    download_cmd = "wget -r -nd -np -A 'viral.*.genomic.fna.gz' ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/"
    # inflate and write to 1 file
    print("Inflating fasta files.")
    extract_cmd = "zcat viral.*.genomic.fna.gz > aux.viruses.fna"

    # changing 'chromosome' names from 'gi|167006425|ref|NC_010314.1|' to NC_010314.1
    sed_cmd = "sed 's/^>gi|[0-9]*|ref|/>/g;s/[;'\'']/_/g;s/|//g' aux.viruses.fna > all.viruses.fna"

    ids_cmd = "grep -oP '[A-Z]{2,2}_[0-9]*\.[0-9]{1,1}' all.viruses.fna > all.viruses.IDs.txt"

    cmds.append(download_cmd)
    cmds.append(extract_cmd)
    cmds.append(sed_cmd)
    cmds.append(ids_cmd)


    q = Queue()
    for cmd in cmds:
        q.submit_nonqueue(cmd, db_path)

def main():
    # build virus reference DB
    virus_db_path = cfg.references["virus_db"]

    IOMethods.create_folder(virus_db_path)

    download_virus_db(virus_db_path)
    
    

if __name__ == "__main__":
    main()
