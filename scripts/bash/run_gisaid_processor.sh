#!/bin/bash

source $1
export COVIGATOR_STORAGE_FOLDER=/scratch/info/projects/covigator/data/gisaid
export COVIGATOR_FORCE_PIPELINE=true

#covigator-processor --source GISAID --local --num-local-cpus 48
covigator-processor --source GISAID --num-jobs 2
