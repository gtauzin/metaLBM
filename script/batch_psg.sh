#!/bin/bash

NPROCS=$1
NTHREADS=1

NTASKS_PER_NODE= $NPROCS

GLOBAL_LENGTH_X=$2
GLOBAL_LENGTH_Y=$2
GLOBAL_LENGTH_Z=$2

LBM_POSTFIX=$3

PARTITION=hsw_p100
NNODES=1
NTASKS_PER_NODE=$(( ${NPROCS}/${NNODES} ))

PARAMS="${NPROCS} ${NTHREADS} ${GLOBAL_LENGTH_X} ${GLOBAL_LENGTH_Y} ${GLOBAL_LENGTH_Z} ${LBM_POSTFIX}"
TARGET_NAME="gpulbm_${NPROCS}_${NTHREADS}_${GLOBAL_LENGTH_X}_${GLOBAL_LENGTH_Y}_${GLOBAL_LENGTH_Z}_${LBM_POSTFIX}"

echo ${PARAMS}
cd ../build
cmake -DPARAMS="${PARAMS}" ..
make ${TARGET_NAME} -j 8

cd ../bin
sbatch --time=4:00:00 --job-name=${TARGET_NAME} \
--partition=${PARTITION} --nodes=${NNODES} --ntasks-per-node=${NTASKS_PER_NODE} \
--mem_bind=local --cpus-per-task=${NTHREADS} \
--output=${TARGET_NAME}-%j.out --error=${TARGET_NAME}-%j.err \
--mail-type=FAIL --mail-user=guillaumetauzin.ut@gmail.com \
--wrap="cd $HOME/Workspace/metaLBM_private/bin && mpirun -np ${NPROCS} ${TARGET_NAME}"
