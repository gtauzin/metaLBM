#!/bin/bash

echo "lbm_stromboli.sh starts with argument: " $1

echo "-- Loading py27"
source activate py27 || true

CMAKEDIR=/home/tauzin/software/cmake
HDF5DIR=/home/tauzin/software/hdf5
FFTWDIR=/home/tauzin/software/fftw
BOOSTDIR=/home/tauzin/software/boost
MPIDIR=/home/tauzin/software/openmpi
GCCDIR=/home/tauzin/software/gcc

export PATH=$PATH:${CMAKEDIR}/bin:${MPIDIR}/bin:${GCCDIR}/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HDF5DIR}/lib:${FFTWDIR}/lib:${BOOSTDIR}/lib:${MPIDIR}/lib:${GCCDIR}/lib:$\
{GCCDIR}/lib64
export CC=${GCCDIR}/bin/gcc
export CXX=${GCCDIR}/bin/g++


#PARTITION=NODE2008
#NODES=8
#NODELIST=node[01,02,03,04,05,06,07,08]
#NODELIST=node[09,10,11,12]
#NTASKSPERNODE=8
#NPROCS_LIST=64
#NTHREADS_LIST=1

PARTITION=NODE2011
NODES=11
#NODELIST=node[13,14,15,16,17,18]
#NODELIST=node[19,20,21,22,23,24]
NODELIST=node[14,15,16,17,18,19,20,21,22,23,24]
NTASKSPERNODE=24
NPROCS_LIST=256
NTHREADS_LIST=1

#PARTITION=NODE2014
#NODES=8
#NODELIST=vulcano[01,02,03,04,05,06,07,08]
#NTASKSPERNODE=32
#NPROCS_LIST=64
#NTHREADS_LIST=1


#BINDING=core:overload-allowed
BINDING=core
MAPPING=core

export OMP_PROC_BIND=true
export OMP_PLACES=cores

DATA_STRUCT=-D_SOA
LOG= #-DENABLE_LOG

echo "-- Initializing simulation on stromboli ${PARTITION} with ${NODES} nodes: ${NODELIST}"
echo "-- Sbatch options: -bind-to ${BINDING}, --map-by mapping ${MAPPING}"
echo "-- OMP options: OMP_PROC_BIND ${OMP_PROC_BIND}, OMP_PLACES ${OMP_PLACES}"

#FILES=$(find ../input/inputJSON -type f -name '*.json')

for NPROCS in $NPROCS_LIST; do
    for NTHREADS in $NTHREADS_LIST; do
        echo "-- NPROCS ${NPROCS}, NTHREADS ${NTHREADS}, DATA_STRUCT ${DATA_STRUCT}"

        echo "-- Relocating to the input folder"
        cd $HOME/Workspace/lbm_solver/input/

        echo "-- Removing input.h and its link"
        rm ../input/inputPy/input.h
        rm ../solver/src/input.h

        echo "-- Removing build files and executable"
        rm -rf ../solver/build/*
        rm ../solver/bin/lbm

        echo "-- Generation input.h from the json file and linking it in src/"
        python inputPy/generate_input_h.py $1
        mv inputPy/input.h ../solver/src/core/

        echo '-- Creating required directories if they don t exist'
        ../script/check_directories.sh solver

        echo "-- Relocating to the build folder to rebuild recompile"
        cd ../solver/build/
        cmake -DMPI=ON -DNPROCS:INT=${NPROCS} -DNTHREADS:INT=${NTHREADS} -DDATA_STRUCT:STRING=${DATA_STRUCT} -DLOG:STRING=${LOG} .. &&
	make lbm -j 8

        echo "-- Relocating to the script folder and adding the job to the job list"
        cd ../../script
        JOBNAME=$(./manage_job.sh lbm $HOSTNAME $1)

        echo "-- Generating the job file and queuing it with sbatch"
        sed -e "s;%PARTITION%;${PARTITION};g" -e "s;%JOBNAME%;${JOBNAME};g" -e "s;%NODES%;${NODES};g" -e "s;%NTASKSPERNODE%;${NTASKSPERNODE};g" -e "s;%NODELIST%;${NODELIST};g" job_stromboli.txt > job_stromboli.sh
        sbatch job_stromboli.sh ${NPROCS} ${BINDING} ${MAPPING}

        echo "-- Removing input.h and its link"
        rm ../input/inputPy/input.h
        rm ../solver/src/input.h

        echo "-- Removing build files and executable"
        rm -rf ../solver/build/*
        #rm ../solver/bin/lbm

    done
done

echo "lbm_stromboli.sh ends"
