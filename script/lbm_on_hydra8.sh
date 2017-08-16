#!/bin/bash

echo "lbm_on_hydra8.sh starts with argument ${1}"

echo "-- Relocating to the script folder"
cd $HOME/Workspace/lbm_solver/script

echo "-- Synchronizing with local repo"
./synchronize.sh hydra8

echo '-- Creating required directories on hydra8 if they don t exist'
sshp hydra8 "cd Workspace/lbm_solver/script/ && ./check_directories.sh"

echo "-- Copying lbm_hydra8.sh file to hydra8"
scpt hydra8 ./lbm_hydra8.sh home/tauzin/Workspace/lbm_solver/script/

echo "-- Copying the json input file to hydra8"
scpt hydra8 $1 home/tauzin/Workspace/lbm_solver/input/inputJSON/

prefix=$(python parse_json.py $1 prefix)
echo "-- Parsing prefix of ${1}: ${prefix}"

echo "-- SSH tp hydra8 to relocate to the script folder there and run run_hydra8.sh"
sshp hydra8 "./Workspace/lbm_solver/script/lbm_hydra8.sh /home/tauzin/Workspace/lbm_solver/input/inputJSON/${prefix}.json"

echo "run_on_hydra8.sh ends"
