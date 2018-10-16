#!/bin/bash

NPROCS=$1
NTHREADS=8

NTASKS_PER_NODE= $NPROCS

GLOBAL_LENGTH_X=$2
GLOBAL_LENGTH_Y=$3
GLOBAL_LENGTH_Z=$4

LBM_POSTFIX=$5

DEP_JOBID=38870

PARAMS="${NPROCS} ${NTHREADS} ${GLOBAL_LENGTH_X} ${GLOBAL_LENGTH_Y} ${GLOBAL_LENGTH_Z} ${LBM_POSTFIX}"
TARGET_NAME="gpulbm_${NPROCS}_${NTHREADS}_${GLOBAL_LENGTH_X}_${GLOBAL_LENGTH_Y}_${GLOBAL_LENGTH_Z}_${LBM_POSTFIX}"

echo ${PARAMS}
cd ../build
cmake -DPARAMS="${PARAMS}" ..
make ${TARGET_NAME} -j 8

cd ../bin

bsub \
    -n ${NPROCS} \
    -gpu "num=4:j_exclusive=yes" \
    -R "span[ptile=4]" \
    -W 24:00 \
    -q normal \
    -oo "${LBM_POSTFIX}.out" \
    -eo "${LBM_POSTFIX}.err" \
    "mpirun -np ${NPROCS} ${TARGET_NAME}"
