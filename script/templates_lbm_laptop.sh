echo "templates_lbm_laptop.sh starts with argument" $1

echo "-- Loading gompi and py27"
module load gompi || true
source activate py27

echo '-- Creating required directories if they don t exist'
cd $HOME/Workspace/lbm_solver/script
./check_directories.sh reference_templates

DATA_STRUCT=-D_SOA

echo "-- Initializing serial simulation on laptop ${PARTITION} with ${NODES} nodes: ${NODELIST}"

#FILES=$(find ../input/inputJSON -type f -name '*.json')

echo "-- DATA_STRUCT ${DATA_STRUCT}"

echo "-- Relocating to the input folder"
cd $HOME/Workspace/lbm_solver/input/

echo "-- Generation input.h from the json file and linking it in solver/src/"
python inputPy/generate_input_h.py $1 'serial'
mv inputPy/input.h ../reference_templates/src/

echo "-- Relocating to the build folder to rebuild and compile"
cd ../reference_templates/build/
cmake ..
make lbm NPROCS=1 NTHREADS=1 DATA_STRUCT=$DATA_STRUCT

echo "-- Relocating to the bin folder and running the job"
cd ../bin/
./lbm

echo "-- Removing input.h and its link"
rm ../../input/inputPy/input.h
rm ../../solver/src/input.h

echo "-- Removing build files and executable"
rm -rf ../build/*
rm lbm

echo "templates_lbm_laptop.sh ends"
