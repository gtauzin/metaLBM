#!/bin/bash

NPROCS=$1
NHOURS=$2
COMMAND=$3

sbatch \
    -n ${NPROCS} \
    -p hsw_p100
    -t ${NHOURS}:00:00 \
    -I \
    --wrap=${COMMAND}
