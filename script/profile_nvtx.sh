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
rm -rf ./CMakeFiles ./CMakeCache.txt

CC=gcc \
CXX=g++ \
CUDACXX=nvcc \
cmake .. \
  -DUSE_CUDA=ON \
  -DUSE_FFTW=ON \
  -DUSE_NVSHMEM=OFF \
  -DUSE_SCOREP=OFF \
  -DUSE_NVTX=ON \
  -DNPROCS=1 \
  -DNTHREADS=1 \
  -DCMAKE_CUDA_FLAGS="-O3 -arch=sm_60" \
  #   -DFFTW_ROOT=<path_to_fftw> \
  -DFFTW_WITH_THREADS=ON \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=On \

cmake -DPARAMS="${PARAMS}" ..
make ${TARGET_NAME} -j 8

cd ../bin
mv ${TARGET_NAME} ${TARGET_NAME}_${POSTFIX}

./interactive_juron 1 " \
LD_PRELOAD=\"/home/gtauzin/software/cuda-profiler/nvtx_pmpi_wrappers/libnvtx_pmpi.so\" mpirun -np ${NPROCS} nvprof -f -o profile_${POSTFIX}.%q{OMPI_COMM_WORLD_RANK}.nvprof ./${TARGET_NAME}_${POSTFIX} \
LD_PRELOAD=\"/home/gtauzin/software/cuda-profiler/nvtx_pmpi_wrappers/libnvtx_pmpi.so\" mpirun -np ${NPROCS} nvprof --analysis-metrics -f -o profile_metrics_${POSTFIX}.%q{OMPI_COMM_WORLD_RANK}.nvprof ./${TARGET_NAME}_${POSTFIX}"
