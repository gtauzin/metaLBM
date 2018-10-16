#!/bin/bash

echo "lbm_on_newturb.sh starts with argument ${1}"

echo "-- Relocating to the script folder"
cd $HOME/Workspace/lbm_solver/script/

echo "-- Synchronizing folder"
./synchronize.sh newturb

echo '-- Creating required directories on stromboli if they don t exist'
sshp newturb "cd /storage/tauzin/Workspace/lbm_solver/script/ && ./check_directories.sh"

echo "-- Copying lbm_newturb.sh file to newturb"
scpt newturb ./lbm_newturb.sh storage/tauzin/Workspace/lbm_solver/script/

echo "-- Copying the json input file to newturb"
scpt newturb $1 storage/tauzin/Workspace/lbm_solver/input/inputJSON/

prefix=$(python parse_json.py $1 prefix)
echo "-- Parsing prefix of ${1}: ${prefix}"

echo "-- SSH to newturb to relocate to the script folder there and run lbm_newturb.sh"
sshp newturb "cd /storage/tauzin/Workspace/lbm_solver/script/ && ./lbm_newturb.sh /storage/tauzin/Workspace/lbm_solver/input/inputJSON/${prefix}.json"

echo "-- SSH to newturb to relocate to the script folder there and run process.sh"
sshp newturb "cd /storage/tauzin/Workspace/lbm_solver/script/ && ./process.sh /storage/tauzin/Workspace/lbm_solver/input/inputJSON/${prefix}.json"

echo "lbm_on_newturb ends.sh"
