echo "lbm_laptop.sh starts with argument" $@

echo "-- Loading gompi and py27"
module load gompi || true
source activate py27

BINDING=core:overload-allowed
MAPPING=core

export OMP_PROC_BIND=true
export OMP_PLACES=cores

BUID_TYPE=Release
NPROCS_LIST=2
NTHREADS_LIST=1
DATA_STRUCT=-D_SOA
LOG= #-DENABLE_LOG

echo "-- Initializing simulation on laptop ${PARTITION} with ${NODES} nodes: ${NODELIST}"
echo "-- MPI options: -bind-to ${BINDING}, --map-by mapping ${MAPPING}"
echo "-- OMP options: OMP_PROC_BIND ${OMP_PROC_BIND}, OMP_PLACES ${OMP_PLACES}"

#FILES=$(find ../input/inputJSON -type f -name '*.json')

echo "-- Looping over listed input JSONFiles"
for JSONFile in "$@"
do
    for NPROCS in $NPROCS_LIST; do
        for NTHREADS in $NTHREADS_LIST; do
            echo "-- NPROCS ${NPROCS}, NTHREADS ${NTHREADS}, DATA_STRUCT ${DATA_STRUCT}"

            echo "-- Relocating to the input folder"
            cd $HOME/Workspace/lbm_solver/input/

            echo "-- Generation input.h from the json file and linking it in solver/src/"
            python inputPy/generate_input_h.py ${JSONFile}
            mv inputPy/input.h ../solver/src/core/

            echo '-- Creating required directories if they don t exist'
            ../script/check_directories.sh solver

            echo "-- Removing previous build files and executable if existing"
            rm -rf ../solver/build/*
            rm ../solver/bin/lbm

            echo "-- Relocating to the build folder to rebuild and compile"
            cd ../solver/build/
            cmake -DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE} -DMPI=ON -DNPROCS:INT=${NPROCS} -DNTHREADS:INT=${NTHREADS} -DDATA_STRUCT:STRING=${DATA_STRUCT} .. &&
                make lbm -j 4

            echo "-- Relocating to the bin folder and running the job"
            cd ../bin/
            mpiexec -np ${NPROCS} ./lbm --bind-to ${BINDING} --map-by ${MAPPING}

            # echo "-- Removing input.h and its link"
            # rm ../../input/inputPy/input.h
            # rm ../src/core/input.h
            # rm ../build/*
            # rm lbm

            # if [[ "$BUILD_TYPE" -eq "Release" ]]
            # then
            #     echo "-- Removing build files and executable"

            # elif [[ "$BUILD_TYPE" -eq "Debug" ]]
            # then
            #     echo "-- Running gprof and removing build files and executable"
            #     gprof ./lbm gmon.out
            #     rm -rf ../build/*
            #     rm lbm
            # fi

        done
    done
done

echo "lbm_laptop.sh ends"
