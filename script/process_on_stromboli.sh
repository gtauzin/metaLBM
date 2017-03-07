#!/bin/bash

echo "--Path of the JSON input file: " $@

echo "-- Synchronizing with local repo"
./synchronize.sh stromboli

echo "-- Looping over listed input JSONFiles"
for JSONFile in "$@"
do
    prefix=$(python parse_json.py $JSONFile prefix)

    scpt stromboli ../input/inputJSON/${prefix}.json /home/tauzin/Workspace/lbm_solver/input/inputJSON/

    sshp stromboli "screen -dm bash -c 'cd Workspace/lbm_solver/script && ./process.sh ../input/inputJSON/${prefix}.json &> process_${prefix}.out'"

done

echo "process.sh ends"
