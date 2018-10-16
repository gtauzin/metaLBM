#!/bin/bash

echo "--Path of the JSON input file: " $1

echo "--Activating python 2.7"
source activate py27

echo "--Parsing the JSON file: "
prefix=$(python parse_json.py $1 prefix)
echo "----prefix: " $prefix
iterationMax=$(python parse_json.py $1 iterationMax)
echo "----iterationMax: " $iterationMax
writeStep=$(python parse_json.py $1 writeStep)
echo "----writeStep: " $writeStep
numberFileMax=$((iterationMax/$writeStep))
echo "----numberFileMax: " $numberFileMax

scp ../input/inputJSON/${prefix}.json num06:/localhdd/Workspace/TurbulentPy/outputPy/

numberFileProcessed=$(ssh num06 "ls -l /home/tauzin/Workspace/TurbulentPy/outputVTR/${prefix}-*" | wc -l)
echo "--Number of files processed: "$numberFileProcessed

while [ "$numberFileProcessed" -lt "$numberFileMax" ]
do
    echo "--There are " $((numberFileMax-$numberFileProcessed)) "files left to process!"
    numberFileLocal=$(ls -l /home/tauzin/Workspace/lbm_solver/output/outputVTR/${prefix}-* | wc -l)
    echo "--Number of files ready: "$numberFileLocal
    sleep 5

    for file in $(seq 0 $numberFileLocal);
    do
        rsync -au --progress "$HOME/Workspace/lbm_solver/output/outputVTR/${prefix}-"${file}.vtr num06:/localhdd/Workspace/TurbulentPy/outputVTR/
    done

    echo "Processing " $((numberFileLocal-$numberFileProcessed)) "files"
    start=$(date +%s.%N)
    ssh num06 "source activate py27; cd /localhdd/Workspace/TurbulentPy/; python process.py /localhdd/Workspace/TurbulentPy/outputPy/${prefix}.json "
    end=$(date +%s.%N)
    diff=$(echo "$end - $start" | bc)
    diff_int=${diff%.*}

    sleep_time=$(echo $(($diff_int>3600?0:3600-$diff_int)) | bc)
    echo "--Sleeping for " $sleep_time "s"
    sleep $sleep_time

    numberFileProcessed=$(ssh num06 "ls -l /home/tauzin/Workspace/TurbulentPy/outputVTR/${prefix}-*" | wc -l)
    echo "--Number of files processed: "$numberFileProcessed
done
