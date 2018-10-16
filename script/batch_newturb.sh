#!/bin/bash

NPROCS=$1
NTHREADS=1

GLOBAL_LENGTH_X=$2
GLOBAL_LENGTH_Y=$2
GLOBAL_LENGTH_Z=$2

LBM_POSTFIX=$3

NNODES=3
NTASKS_PER_NODE=24
#NTASKS_PER_NODE=$(( ${NPROCS}/${NNODES} ))

PARAMS="${NPROCS} ${NTHREADS} ${GLOBAL_LENGTH_X} ${GLOBAL_LENGTH_Y} ${GLOBAL_LENGTH_Z} ${LBM_POSTFIX}"
TARGET_NAME="cpulbm_${NPROCS}_${NTHREADS}_${GLOBAL_LENGTH_X}_${GLOBAL_LENGTH_Y}_${GLOBAL_LENGTH_Z}_${LBM_POSTFIX}"

echo ${PARAMS}
cd ../build
cmake -DPARAMS="${PARAMS}" ..
make VERBOSE=1 ${TARGET_NAME} -j 8


# cd ../bin
# echo "cd /mnt/beegfs/tauzin/metaLBM_private/bin && mpirun -np ${NPROCS} ./${TARGET_NAME}" \
# | qsub -N ${TARGET_NAME} \
# -l nodes=${NNODES}:ppn=${NTASKS_PER_NODE},walltime=24:00:00 \
# -o stdout_${TARGET_NAME} -e stderr_${TARGET_NAME} \
# -m abe -M guillaumetauzin.ut@gmail.com

#echo "NOT SUBMITTED"
