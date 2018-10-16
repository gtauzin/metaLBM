#!/bin/bash

NPROCS=$1
NHOURS=$2

bsub \
    -n ${NPROCS} \
    -gpu "num=4:j_exclusive=yes" \
    -R "span[ptile=4]" \
    -W ${NHOURS}:00 \
    -x \
    -Is \
    -tty /bin/bash
