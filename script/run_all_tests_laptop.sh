echo "run_all_tests_laptop.sh starts"

echo "-- Loading gompi and py27"
module load gompi || true
source activate py27

echo "-- Relocating to the script folder"
cd $HOME/Workspace/lbm_solver/script/

BINDING=core:overload-allowed
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

echo "-- Relocating to the input folder"
cd $HOME/Workspace/lbm_solver/input/

echo "-- Generation input.h from the json file and linking it in solver/src/"
python inputPy/generate_input_h.py $1
mv inputPy/input.h ../solver/src/core/

echo '-- Creating required directories if they don t exist'
../script/check_directories.sh solver

echo "-- Removing previous build files and executable if existing"
rm -rf ../solver/build/*
rm ../solver/bin/all_tests

echo "-- Relocating to the build folder to rebuild and compile"
cd ../solver/build/
cmake -DMPI=ON -DNPROCS:INT=${NPROCS} -DNTHREADS:INT=${NTHREADS} -DDATA_STRUCT:STRING=${DATA_STRUCT} .. &&
make all_tests

echo "-- Relocating to the bin folder and running the job"
cd ../bin/
mpiexec -np ${NPROCS} ./all_tests --log_level=test_suite

echo "-- Removing input.h and its link"
rm ../../input/inputPy/input.h
rm ../src/core/input.h
rm -rf ../build/*
rm all_tests

echo "run_all_tests_laptop.sh ends"
