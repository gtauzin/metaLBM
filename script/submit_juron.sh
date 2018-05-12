#!/bin/bash

NPROCS=$1
NTHREADS=1

NTASKS_PER_NODE= $NPROCS

GLOBAL_LENGTH_X=$2
GLOBAL_LENGTH_Y=$2
GLOBAL_LENGTH_Z=$2

POSTFIX=$3

PARAMS="${NPROCS} ${NTHREADS} ${GLOBAL_LENGTH_X} ${GLOBAL_LENGTH_Y} ${GLOBAL_LENGTH_Z}"
TARGET_NAME="gpulbm_${NPROCS}_${NTHREADS}_${GLOBAL_LENGTH_X}_${GLOBAL_LENGTH_Y}_${GLOBAL_LENGTH_Z}"

echo ${PARAMS}
cd ../build
cmake -DPARAMS="${PARAMS}" ..
make ${TARGET_NAME} -j 8

cd ../bin
mv ${TARGET_NAME} ${POSTFIX}_${TARGET_NAME}

bsub \
    -n ${NPROCS} \
    -x \
    -R "span[ptile=${NPROCS}]" \
    -R "rusage[ngpus_shared=${NPROCS}]" \
    -W 24:00 \
    -q normal \
    -o ${POSTFIX}.out \
    -e ${POSTFIX}.err \
    "mpirun -np ${NPROCS} ${POSTFIX}_${TARGET_NAME}"
