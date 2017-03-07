#!/bin/bash

echo "plot.sh starts with arguments: " $@

echo "-- Loading python 2.7"
source activate py27

echo "-- Looping over listed input JSONFiles"
for JSONFile in "$@"
do
    echo "-- Relocating to the script folder"
    cd ../script
    echo "-- Parsing the JSON file: " $JSONFile
    prefix=$(python parse_json.py $JSONFile prefix)

    echo "-- Plotting "
    cd ../processing/
    python plot.py ../input/inputJSON/${prefix}.json
done

echo "plot.sh ends"
