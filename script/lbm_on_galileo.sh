#!/bin/bash

echo "lbm_on_galileo.sh starts with argument ${1}"

echo "-- Relocating to the script folder"
cd $HOME/Workspace/lbm_solver/script/

echo "-- Synchronizing folder"
./synchronize.sh galileo

echo '-- Creating required directories on stromboli if they don t exist'
sshp galileo "cd /storage/tauzin/Workspace/lbm_solver/script/ && ./check_directories.sh"

echo "-- Copying lbm_galileo.sh file to galileo"
scpt galileo ./lbm_galileo.sh gpfs/scratch/userexternal/gtauzin0/Workspace/lbm_solver/script/

echo "-- Copying the json input file to galileo"
scpt galileo $1 gpfs/scratch/userexternal/gtauzin0/Workspace/lbm_solver/input/inputJSON/

prefix=$(python parse_json.py $1 prefix)
echo "-- Parsing prefix of ${1}: ${prefix}"

echo "-- SSH to galileo to relocate to the script folder there and run lbm_galileo.sh"
sshp galileo "cd /gpfs/scratch/userexternal/gtauzin0/Workspace/lbm_solver/script/ && ./lbm_galileo.sh /gpfs/scratch/userexternal/gtauzin0/Workspace/lbm_solver/input/inputJSON/${prefix}.json"

echo "lbm_on_galileo ends.sh"
