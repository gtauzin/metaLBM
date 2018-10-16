#!/bin/bash

echo "lbm_on_stromboli.sh starts with argument ${1}"

echo "-- Relocating to the script folder"
cd $HOME/Workspace/lbm_solver/script/

echo "-- Synchronizing with local repo"
./synchronize.sh stromboli

echo '-- Creating required directories on stromboli if they don t exist'
sshp stromboli "cd Workspace/lbm_solver/script/ && ./check_directories.sh"

echo "-- Copying lbm_stromboli.sh file to stromboli"
scpt stromboli ./lbm_stromboli.sh home/tauzin/Workspace/lbm_solver/script/

echo "-- Copying the json input file to stromboli"
scpt stromboli $1 home/tauzin/Workspace/lbm_solver/input/inputJSON/

prefix=$(python parse_json.py $1 prefix)
echo "-- Parsing prefix of ${1}: ${prefix}"

echo "-- SSH tp stromboli to relocate to the script folder there and run lbm_stromboli.sh"
sshp stromboli "cd Workspace/lbm_solver/script/ && ./lbm_stromboli.sh /home/tauzin/Workspace/lbm_solver/input/inputJSON/${prefix}.json"

echo "-- Running process_on_stromboli.sh"
./process.sh ../input/inputJSON/${prefix}.json"

echo "lbm_on_stromboli.sh ends"
