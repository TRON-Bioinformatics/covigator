#!/usr/bin/env python

import getpass
import os

pipeline_name = "COVIGATOR"
tools = ("BWA", "PileupVCF", "Collect")
sender = "CoVigator@TrOn-Mainz.DE"
receiver = ("Patrick.Sorn@TrOn-Mainz.DE")
queueing_system = "slurm"
time_limit = "30-00:00:0"
user = getpass.getuser()
partition = "Compute"
module_dir = os.path.dirname(os.path.realpath(__file__))
#fastadir = "/code/plasmidify/data/fasta/"


paths = {
    "novo_folder": "/scratch/info/projects/SARS-CoV-2/Novoalign/results",
    "uncompressed_fastqs": "/scratch/info/projects/SARS-CoV-2/Novoalign/fastqs",
}


cmds = {
    "virgene": "/scratch/info/projects/CM23_VirusID/VIRGENE_1.0/VIRGENE/VIRGENE.py",
    "samtools": "/code/samtools/1.9/samtools",
    "msa": os.path.join(module_dir, "surface_glyco.r"),
    "bwa": "/code/bwa/0.7.17/bwa",
    "bcftools": "/code/bcftools/1.9/bcftools",
    "snpeff": "java -jar /code/snpEFF_latest_core/snpEff/snpEff.jar",
}



references = {
    #"bwa_ref": "/",
    "sra_fastq_folder": "/scratch/info/projects/SARS-CoV-2/sra",
    "sra2_fastq_folder": "/scratch/info/projects/SARS-CoV-2/sra2",
    "novo_index": "/scratch/info/projects/SARS-CoV-2/Novoalign/novoindex/MN908947.3.nix",
    "novo_fasta": "/scratch/info/projects/SARS-CoV-2/Novoalign/novoindex/MN908947.3.fa",
    "novo_bed": "/scratch/info/projects/SARS-CoV-2/Novoalign/novoindex/MN908947.3_gp02_Sgene.bed",
    "virus_db": os.path.join(module_dir, "viruses_index"),
    #"virus_index": "/projects/data/virus/viral_db_20200401/viruses_index",
}

