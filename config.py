
import os

paths = {
    "novo_folder": "/scratch/info/projects/SARS-CoV-2/Novoalign/results",
    "uncompressed_fastqs": "/scratch/info/projects/SARS-CoV-2/Novoalign/fastqs",
}


cmds = {
    "virgene": "/scratch/info/projects/CM23_VirusID/VIRGENE_1.0/VIRGENE/VIRGENE.py",
}

module_dir = os.path.dirname(os.path.realpath(__file__))

references = {
    "sra_fastq_folder": "/scratch/info/projects/SARS-CoV-2/sra",
    "sra2_fastq_folder": "/scratch/info/projects/SARS-CoV-2/sra2",
    "novo_index": "/scratch/info/projects/SARS-CoV-2/Novoalign/novoindex/MN908947.3.nix",
    "novo_fasta": "/scratch/info/projects/SARS-CoV-2/Novoalign/novoindex/MN908947.3.fa",
    "novo_bed": "/scratch/info/projects/SARS-CoV-2/Novoalign/novoindex/MN908947.3_gp02_Sgene.bed",
    "virus_db": os.path.join(module_dir, "viruses_index"),
    #"virus_index": "/projects/data/virus/viral_db_20200401/viruses_index",
}

