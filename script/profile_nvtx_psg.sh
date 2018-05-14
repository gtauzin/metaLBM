#!/bin/bash

NPROCS=$1
NTHREADS=1

GLOBAL_LENGTH_X=$2
GLOBAL_LENGTH_Y=$2
GLOBAL_LENGTH_Z=$2

POSTFIX=$3

PARAMS="${NPROCS} ${NTHREADS} ${GLOBAL_LENGTH_X} ${GLOBAL_LENGTH_Y} ${GLOBAL_LENGTH_Z}"
TARGET_NAME="gpulbm_${NPROCS}_${NTHREADS}_${GLOBAL_LENGTH_X}_${GLOBAL_LENGTH_Y}_${GLOBAL_LENGTH_Z}"

cd ../build
./build.sh
cmake -DUSE_NVTX=ON -DPARAMS="${PARAMS}" ..
make ${TARGET_NAME} -j 8

cd ../bin
mv ${TARGET_NAME} ${TARGET_NAME}_${POSTFIX}

cd ../script
./interactive_psg.sh ${NPROCS} 1 " \
LD_PRELOAD=\"${NVTX_MPI_WRAPPERS}/libnvtx_pmpi.so\" mpirun -np ${NPROCS} nvprof -f -o profile_${POSTFIX}.%q{OMPI_COMM_WORLD_RANK}.nvprof ../bin/${TARGET_NAME}_${POSTFIX}; \ 
LD_PRELOAD=\"${NVTX_MPI_WRAPPERS}/libnvtx_pmpi.so\" mpirun -np ${NPROCS} nvprof --analysis-metrics -f -o profile_metrics_${POSTFIX}.%q{OMPI_COMM_WORLD_RANK}.nvprof ../bin/${TARGET_NAME}_${POSTFIX}"
