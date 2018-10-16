#!/bin/bash                                                                                                           

mkdir ../build
cd ../build 
rm -rf ./CMakeFiles ./CMakeCache.txt CTestTestfile.cmake Makefile cmake_install.cmake

CC=gcc \
CXX=g++ \
CUDACXX=nvcc \
cmake .. \
   -DUSE_CUDA=ON \
   -DUSE_FFTW=ON \
   -DUSE_NVSHMEM=OFF \
   -DUSE_SCOREP=OFF \
   -DUSE_NVTX=OFF \
   -DNPROCS=1 \
   -DNTHREADS=1 \
   -DCMAKE_CUDA_FLAGS="-O3 -arch=sm_30" \
#"-O3 -arch=sm_30" \                                                                                                  
#   -DFFTW_ROOT=<path_to_fftw> \                                                                                      
#   -DFFTW_WITH_THREADS=ON \                                                                                          
#   -DCMAKE_EXPORT_COMPILE_COMMANDS=On \                                                                              

# Debug flags gpu : -O0 -lineinfo -G                                                                                  
# Debug flags cpu : -O2 -g -pg     
