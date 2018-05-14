#!/bin/bash

NPROCS=$1
TARGET_NAME="gpulbm_${NPROCS}_${NTHREADS}_${GLOBAL_LENGTH_X}_${GLOBAL_LENGTH_Y}_${GLOBAL_LENGTH_Z}"

mpirun -np ${NPROCS} tau_exec -T mpi,cupti -ebs ../bin/${TARGET_NAME}
