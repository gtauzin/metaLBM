#!/bin/bash
. ~/.bashrc

echo "fetch_and_plot_num06.sh starts with argument: " $1

echo "-- Parsing the JSON file: "
prefix=$(python parse_json.py $1 prefix)
echo "---- prefix: " $prefix

echo "-- Fetching files from num06"
scp num06:/localhdd/Workspace/lbm_solver/output/outputVTR/${prefix}_processed_VTRfiles.txt ../output/outputPy/
scp num06:/localhdd/Workspace/lbm_solver/output/outputPy/${prefix}* ../output/outputPy/

echo "-- Loading python 2.7"
source activate py27

echo "-- Plotting"
cd ../processing
python plot.py ../output/outputPy/${prefix}.json

echo "fetch_and_plot_num06.sh ends"
