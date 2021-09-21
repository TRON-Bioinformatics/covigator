#!/bin/bash

source $1
export COVIGATOR_FORCE_PIPELINE=true
# /projects/SARS-CoV-2/gisaid/sequences_fasta_2021_08_17/sequences.fasta
fasta=$2
# /projects/SARS-CoV-2/gisaid/metadata_tsv_2021_08_17/metadata.tsv
metadata=$3

covigator-gisaid-accessor --input-fasta $fasta --input-metadata $metadata
