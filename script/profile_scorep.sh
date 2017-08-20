#!/bin/bash

cd ~/Workspace/metaLBM_private/
mkdir build
cd build
rm -rf *
SCOREP_WRAPPER=OFF CC=scorep-gcc CXX=scorep-g++ cmake .. -DCMAKE_BUILD_TYPE=Debug
SCOREP_WRAPPER=ON make lbm
cd src

export SCOREP_ENABLE_TRACING=true
mpirun -np $1 ./lbm
