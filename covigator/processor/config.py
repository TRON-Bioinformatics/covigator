#!/usr/bin/env python

import getpass
import os

pipeline_name = "COVIGATOR"
tools = ("BWA", "LoFreq", "SNPeff")

module_dir = os.path.dirname(os.path.realpath(__file__))

resources = {
    "bwa": {
        "cpu": 6,
        "mem": 10
    },
    "pileupvcf": {
        "cpu": 1,
        "mem": 10
    },
    "lofreq": {
        "cpu": 1,
        "mem": 10
    },
    "snpeff": {
        "cpu": 1,
        "mem": 10
    }
}

# define path to input and output folders
paths = {
    "novo_folder": "/scratch/info/projects/SARS-CoV-2/Novoalign/results",
    "uncompressed_fastqs": "/scratch/info/projects/SARS-CoV-2/Novoalign/fastqs",
}

# paths to executables
cmds = {
    "virgene": "/scratch/info/projects/CM23_VirusID/VIRGENE_1.0/VIRGENE/VIRGENE.py",
    "samtools": "/code/samtools/1.9/samtools",
    "msa": os.path.join(module_dir, "surface_glyco.r"),
    "bwa": "/code/bwa/0.7.17/bwa",
    "bcftools": "/code/bcftools/1.9/bcftools",
    "lofreq": "/code/lofreq_star-2.0.0-beta-3/lofreq/lofreq",
    "snpeff": "java -jar /code/snpEFF_latest_core/snpEff/snpEff.jar",
}


# paths to reference files
references = {
    "sra_fastq_folder": "/scratch/info/projects/SARS-CoV-2/sra",
    "sra2_fastq_folder": "/scratch/info/projects/SARS-CoV-2/sra2",
    "novo_index": "/scratch/info/projects/SARS-CoV-2/index/MN908947.3.nix",
    "novo_fasta": "/scratch/info/projects/SARS-CoV-2/index/MN908947.3.fa",
    "novo_bed": "/scratch/info/projects/SARS-CoV-2/Novoalign/novoindex/MN908947.3_gp02_Sgene.bed",
    "virus_db": os.path.join(module_dir, "viruses_index"),
}

