#!/bin/bash

NHOURS=$1
COMMAND=$2

bsub \
    -n 1 \
    -R "span[ptile=4]" \
    -R "rusage[ngpus_shared=4]" \
    -W ${NHOURS}:00 \
    -x \
    -Is \
    -E "${COMMAND}"
