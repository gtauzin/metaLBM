#!/bin/bash

NPROCS=$1
NHOURS=$2

bsub \
    -n ${NPROCS} \
    -R "span[ptile=4]" \
    -R "rusage[ngpus_shared=4]" \
    -W ${NHOURS}:00 \
    -x \
    -Is \
    -tty /bin/bash
