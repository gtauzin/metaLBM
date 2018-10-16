#!/bin/bash

echo "lbm_hydra8.sh starts with argument: " $1

echo "-- Loading gompi and py27"
source activate py27 || true

BOOSTDIR=/home/tauzin/software
MPIDIR=/usr/local/openmpi
CUDADIR=/usr/local/cuda-7.5

export CUDA_INCLUDE=${CUDADIR}/include
export CUDA_LIB64=${CUDADIR}/lib64

export MPI_BIN=${MPIROOT}/bin
export MPI_INCLUDE=${MPIROOT}/include
export MPI_LIB=${MPIROOT}/lib

#export LD_LIBRARY_PATH=${CUDA_LIB64}:${MPI_LIB}:${LD_LIBRARY_PATH}
#export PATH=${MPI_BIN}:${PATH}

export PATH=/usr/local/cuda-7.5/bin:/usr/local/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda-7.5/lib64:/usr/local/openmpi/lib:/usr/lib:/usr/local/openmpi/lib64:/home/tauzin/software/lib:$LD_LIBRARY_PATH
#export PATH=$PATH:${MPIDIR}/bin:${GCCDIR}/bin
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${BOOSTDIR}/lib:${MPIDIR}/lib:${GCCDIR}/lib:${GCCDIR}/lib64

#BINDING=core:overload-allowed
BINDING=core
MAPPING=core

export OMP_PROC_BIND=true
export OMP_PLACES=cores

NPROCS_LIST=2
NTHREADS_LIST=1
DATA_STRUCT=-D_SOA

echo "-- Initializing simulation on hydra8 ${PARTITION} with ${NODES} nodes: ${NODELIST}"
echo "-- Sbatch options: -bind-to ${BINDING}, --map-by mapping ${MAPPING}"
echo "-- OMP options: OMP_PROC_BIND ${OMP_PROC_BIND}, OMP_PLACES ${OMP_PLACES}"

#FILES=$(find ../input/inputJSON -type f -name '*.json')

for NPROCS in $NPROCS_LIST; do
    for NTHREADS in $NTHREADS_LIST; do
        echo "-- NPROCS ${NPROCS}, NTHREADS ${NTHREADS}, DATA_STRUCT ${DATA_STRUCT}"

        echo "-- Relocating to the input folder"
        cd $HOME/Workspace/lbm_solver/input/

        echo "-- Generation input.h from the json file and linking it in src/"
        python inputPy/generate_input_h.py $1
        mv inputPy/input.h ../solver/src/core/

        echo '-- Creating required directories if they don t exist'
        ../script/check_directories.sh solver

        echo "-- Removing previous build files and executable"
        rm -rf ../solver/build/*
        rm ../solver/bin/lbm

        echo "-- Relocating to the build folder to rebuild recompile"
        cd ../solver/build/
        cmake -DMPI=ON -DNPROCS:INT=${NPROCS} -DNTHREADS:INT=${NTHREADS} -DDATA_STRUCT:STRING=${DATA_STRUCT} .. &&
	make lbm

        echo "-- Relocating to the script folder and adding the job to the job list"
        cd ../../script
        JOBNAME=$(./manage_job.sh lbm $HOSTNAME $1)

        echo "-- Generating the job file and queuing it with sbatch"
        cd ../solver/bin
        #mpiexec -np ${NPROCS} --map-by ${BINDING} -bind-to ${MAPPING} ./lbm
        mpiexec -np ${NPROCS} ./lbm

        echo "-- Removing input.h and its link"
        rm ../input/inputPy/input.h
        rm ../solver/src/core/input.h

        echo "-- Removing build files and executable"
        rm -rf ../solver/build/*
        rm ../solver/bin/lbm

    done
done

echo "lbm_hydra8.sh ends"
