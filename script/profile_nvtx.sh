#!/bin/bash

NPROCS=$1
TARGET_NAME=$2

LD_PRELOAD="${NVTX_MPI_WRAPPERS}/libnvtx_pmpi.so" mpirun -np ${NPROCS} \
  nvprof -f -o profile_${TARGET_NAME}.%q{OMPI_COMM_WORLD_RANK}.nvprof \
  ../bin/${TARGET_NAME}  

LD_PRELOAD="${NVTX_MPI_WRAPPERS}/libnvtx_pmpi.so" mpirun -np ${NPROCS} \
  nvprof --analysis-metrics -f -o profile_metrics_${TARGET_NAME}.%q{OMPI_COMM_WORLD_RANK}.nvprof \
  ../bin/${TARGET_NAME}
