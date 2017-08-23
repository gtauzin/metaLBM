#!/bin/bash

cd ~/metaLBM_private/
rm -rf build
mkdir build
cd build
SCOREP_WRAPPER=OFF cmake .. -DCMAKE_C_COMPILER=scorep-mpicc -DCMAKE_CXX_COMPILER=scorep-mpic++
export SCOREP_WRAPPER=ON
export SCOREP_ENABLE_TRACING=true
export SCOREP_TOTAL_MEMORY=3g

make lbm VERBOSE=1 SCOREP_WRAPPER_INSTRUMENTER_FLAGS="--verbose --user --nocompiler"

cd ../script
./submit_job.sh
