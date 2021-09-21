#!/bin/bash

source $1

covigator-processor --num-jobs 1 --source ENA
#covigator-processor --source ENA --local --num-local-cpus 5
