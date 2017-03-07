#!/bin/bash

echo "run_all_tests_stromboli.sh starts"

echo "-- Loading gompi and py27"
module load gompi || true
module load CMake || true
source activate py27

BOOSTDIR=/home/tauzin/software
MPIDIR=/home/tauzin/software/easybuild/software/OpenMPI/1.10.2-GCC-4.9.3-2.25
GCCDIR=/home/tauzin/software/easybuild/software/GCCcore/4.9.3

export PATH=$PATH:${MPIDIR}/bin:${GCCDIR}/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${BOOSTDIR}/lib:${MPIDIR}/lib:${GCCDIR}/lib:${GCCDIR}/lib64
export CC=${GCCDIR}/bin/gcc
export CXX=${GCCDIR}/bin/g++

BINDING=core
MAPPING=core

export OMP_PROC_BIND=true
export OMP_PLACES=cores

NPROCS=1
NTHREADS=1
DATA_STRUCT=-D_SOA

echo "-- Initializing simulation on laptop ${PARTITION} with ${NODES} nodes: ${NODELIST}"
echo "-- MPI options: -bind-to ${BINDING}, --map-by mapping ${MAPPING}"
echo "-- OMP options: OMP_PROC_BIND ${OMP_PROC_BIND}, OMP_PLACES ${OMP_PLACES}"

echo "-- NPROCS ${NPROCS}, NTHREADS ${NTHREADS}, DATA_STRUCT ${DATA_STRUCT}"

echo "-- Relocating to the script folder"
cd $HOME/Workspace/lbm_solver/script/

echo "-- Linking test input.h it in src/"
cp ../solver/test/input.h ../solver/src/

echo "-- Removing previously generated executable files"
rm ../solver/bin/all_tests > /dev/null

echo '-- Creating required directories if they don t exist'
./check_directories.sh

echo "-- Relocating to the build folder to rebuild and compile"
cd ../solver/build/
cmake .. > /dev/null
make all_tests NPROCS=${NPROCS} NTHREADS=${NTHREADS} DATA_STRUCT=$DATA_STRUCT > /dev/null

echo "-- Relocating to the bin folder and running the job"
cd ../bin/
mpiexec -np ${NPROCS} all_tests --log_level=test_suite

echo "-- Removing input.h from src/"
rm ../src/input.h

echo "-- Removing build files and executables"
rm -rf ../build/* > /dev/null

echo "run_all_tests_stromboli.sh ends"
