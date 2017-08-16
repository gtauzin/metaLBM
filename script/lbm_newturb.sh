#!/bin/bash

echo "lbm_newturb.sh starts with argument: " $1

echo "-- Loading py27"
#export PATH="/scratch/tauzin/software/anaconda3/envs/py27/bin:$PATH"
source /scratch/tauzin/software/anaconda3/envs/py27/bin/activate py27

BOOSTDIR=/scratch/tauzin/software
VTKDIR=/scratch/tauzin/software/VTK-build
MPIDIR=/opt/openmpi/1.10.1

export PATH=$PATH:${MPIDIR}/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MPIDIR}/lib:${BOOSTDIR}/lib:${VTKDIR}/lib
export CC=${GCCDIR}/bin/gcc
export CXX=${GCCDIR}/bin/g++

#BINDING=core:overload-allowed
BINDING=core
MAPPING=core

export OMP_PROC_BIND=true
export OMP_PLACES=cores

#NODES=11
#NTASKSPERNODE=24
#NPROCS_LIST=256
#NTHREADS_LIST=1

NODES=6
NTASKSPERNODE=24
NPROCS_LIST=128
NTHREADS_LIST=1


#NODES=22
#NTASKSPERNODE=24
#NPROCS_LIST=512
#NTHREADS_LIST=1

#NODES=1
#NTASKSPERNODE=1
#NPROCS_LIST=1
#NTHREADS_LIST=1


DATA_STRUCT=-D_SOA

echo "-- Initializing simulation on newturb ${PARTITION} with ${NODES} nodes: ${NODELIST}"
echo "-- Sbatch options: -bind-to ${BINDING}, --map-by mapping ${MAPPING}"
echo "-- OMP options: OMP_PROC_BIND ${OMP_PROC_BIND}, OMP_PLACES ${OMP_PLACES}"

#FILES=$(find ../input/inputJSON -type f -name '*.json')

for NPROCS in $NPROCS_LIST; do
    for NTHREADS in $NTHREADS_LIST; do
        echo "-- NPROCS ${NPROCS}, NTHREADS ${NTHREADS}, DATA_STRUCT ${DATA_STRUCT}"

        echo "-- Relocating to the input folder"
        cd /storage/tauzin/Workspace/lbm_solver/input/

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
        cmake -DMPI_OMP=ON -DNPROCS:INT=${NPROCS} -DNTHREADS:INT=${NTHREADS} -DDATA_STRUCT:STRING=${DATA_STRUCT} .. &&
	make lbm

        echo "-- Relocating to the script folder and adding the job to the job list"
        cd ../../script
        JOBNAME=$(./manage_job.sh lbm $HOSTNAME $1)

        echo "-- Generating the job file and queuing it with sbatch"
        sed -e "s;%JOBNAME%;${JOBNAME};g" -e "s;%NODES%;${NODES};g" -e "s;%NPROCS%;${NPROCS};g" -e "s;%NTASKSPERNODE%;${NTASKSPERNODE};g" job_newturb.txt > job_newturb.psb
        qsub job_newturb.psb #${BINDING} ${MAPPING}

        echo "-- Removing input.h and its link"
        rm ../input/inputPy/input.h
        rm ../solver/src/input.h

        echo "-- Removing build files"
        rm -rf ../solver/build/*

    done
done

echo "lbm_newturb.sh ends"
