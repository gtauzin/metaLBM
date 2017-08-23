#!/bin/bash

NPROCS=2
NTHREADS=1
NODES=1
NTASKSPERNODE=2

JOBNAME="lbm_${NPROCS}_${NTHREADS}"
echo "-- Generating the job file and queuing it with sbatch"
sed -e "s;%JOBNAME%;${JOBNAME};g" -e "s;%NODES%;${NODES};g" -e "s;%NTASKSPERNODE%;${NTASKSPERNODE};g" -e "s;%NPROCS%;${NPROCS};g" job.txt > job.sh
sbatch job.sh
