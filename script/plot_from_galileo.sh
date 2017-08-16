#!/bin/bash

echo "plot_from_galileo.sh starts with argument ${1}"

echo "-- Loading python 2.7"
source activate py27

echo "-- Relocating to the script folder"
cd $HOME/Workspace/lbm_solver/script/

prefix=$(python parse_json.py $1 prefix)
echo "-- Parsing prefix of ${1}: ${prefix}"

echo "-- Copying processed files to local directory"
scpf galileo /gpfs/scratch/userexternal/gtauzin0/Workspace/lbm_solver/output/outputPy/${prefix}* ../output/outputPy/

echo "-- Relocating to the processing folder and run plot.py"
cd ../processing && python plot.py ../input/inputJSON/${prefix}.json

echo "-- Opening generated plots in nautilus"
nautilus ../output/outputPlot

echo "lbm_on_galileo.sh ends"
