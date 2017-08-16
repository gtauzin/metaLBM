#!/bin/bash

echo "lbm_galileo.sh starts with argument: " $1

echo "-- Loading modules and py27"
source /galileo/home/userexternal/gtauzin0/Software/anaconda2/bin/activate py27
module load cmake

module load autoload boost/1.58.0--openmpi--1.8.4--gnu--4.9.2

#export MPICC=mpicc
#export OMPI_CC=icc

BOOSTDIR=/galileo/home/userexternal/gtauzin0/software
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${BOOSTDIR}/lib

#BINDING=core:overload-allowed
BINDING=socket
MAPPING=socket

PROJECT=INF17_fldturb

NODES=8
NCPUS=16
MPIPROCSPERNODE=16
NPROCS_LIST=128
NTHREADS_LIST=1

#NODES=1
#NCPUS=1
#MPIPROCSPERNODE=1
#NPROCS_LIST=1
#NTHREADS_LIST=1

#NODES=4
#NCPUS=16
#MPIPROCSPERNODE=4
#NPROCS_LIST=16
#NTHREADS_LIST=1

MEM=118GB
#WALLTIME=00:00:30

DATA_STRUCT=-D_SOA
LOG= #-DENABLE_LOG

echo "-- Initializing simulation on galileo ${PARTITION} with ${NODES} nodes: ${NODELIST}"
echo "-- Sbatch options: -bind-to ${BINDING}, --map-by mapping ${MAPPING}"
echo "-- OMP options: OMP_PROC_BIND ${OMP_PROC_BIND}, OMP_PLACES ${OMP_PLACES}"

#FILES=$(find ../input/inputJSON -type f -name '*.json')

for NPROCS in $NPROCS_LIST; do
    for NTHREADS in $NTHREADS_LIST; do
        echo "-- NPROCS ${NPROCS}, NTHREADS ${NTHREADS}, DATA_STRUCT ${DATA_STRUCT}"

        echo "-- Relocating to the input folder"
        cd /gpfs/scratch/userexternal/gtauzin0/Workspace/lbm_solver/input/

        echo "-- Generation input.h from the json file and linking it in src/"
        python inputPy/generate_input_h.py $1
        mv inputPy/input.h ../solver/src/core/

        echo '-- Creating required directories if they don t exist'
        ../script/check_directories.sh solver

        echo "-- Removing previous executable and build files if existing"
        rm ../solver/bin/lbm
        rm -rf ../solver/build/*

        echo "-- Relocating to the build folder to rebuild recompile"
        cd ../solver/build/
        cmake -DICC=OFF -DMPI=ON -DNPROCS:INT=${NPROCS} -DNTHREADS:INT=${NTHREADS} -DDATA_STRUCT:STRING=${DATA_STRUCT} -DLOG:STRING=${LOG} .. &&
	make lbm -j 8

        echo "-- Relocating to the script folder and adding the job to the job list"
        cd ../../script
        JOBNAME=$(./manage_job.sh lbm $HOSTNAME $1)

        echo "-- Generating the job file and queuing it with sbatch"
        sed -e "s;%PROJECT%;${PROJECT};g" -e "s;%JOBNAME%;${JOBNAME};g" -e "s;%NODES%;${NODES};g" -e "s;%NPROCS%;${NPROCS};g" -e "s;%MAPPING%;${MAPPING};g" -e "s;%BINDING%;${BINDING};g" -e "s;%MPIPROCSPERNODE%;${MPIPROCSPERNODE};g" -e "s;%NCPUS%;${NCPUS};g" -e "s;%MEM%;${MEM};g" -e "s;%WALLTIME%;${WALLTIME};g" job_galileo.txt > job_galileo.psb
        qsub job_galileo.psb #${BINDING} ${MAPPING}

        echo "-- Removing input.h and its link"
        rm ../input/inputPy/input.h
        rm ../solver/src/input.h

        echo "-- Removing build files"
        rm -rf ../solver/build/*

    done
done

echo "lbm_galileo.sh ends"
