#!/bin/bash

echo "process.sh starts with arguments: " $@

echo "-- Loading python 2.7"
source activate py27

echo "-- Looping over listed input JSONFiles"
for JSONFile in "$@"
do
    echo "-- Relocating to the script folder"
    cd ../script
    echo "-- Parsing the JSON file: " $JSONFile
    prefix=$(python parse_json.py $JSONFile prefix)
    echo "---- prefix: " $prefix
    iterationMax=$(python parse_json.py $JSONFile iterationMax)
    echo "---- iterationMax: " $iterationMax
    writeStep=$(python parse_json.py $JSONFile writeStep)
    echo "---- writeStep: " $writeStep
    numberFileMax=$((iterationMax/$writeStep))-2
    echo "---- numberFileMax: " $numberFileMax

    numberFileProcessed=0
    echo "-- Number of files processed: "$numberFileProcessed

    while [[ "$numberFileProcessed" -lt "$numberFileMax" ]]
    do
        echo "-- There are " $((numberFileMax-$numberFileProcessed)) "files left to process!"
        numberFileLocal=$(ls -l ../output/outputVTR/${prefix}-* | wc -l)
        echo "-- Number of files ready: "$numberFileLocal
        sleep 5

        echo "-- Processing " $((numberFileLocal-$numberFileProcessed)) "files"
        cd ../processing/
        python process.py ../input/inputJSON/${prefix}.json

        numberFileProcessed=$numberFileLocal
        echo "-- Number of files processed: "$numberFileProcessed
        if [[ "$numberFileProcessed" -gt "$numberFileMax" ]]
        then
            echo "-- All files have been processed!"
            continue
        else
            sleep_time=180
            echo "-- Sleeping for " $sleep_time "s"
            sleep $sleep_time
        fi

    done
done
echo "process.sh ends"
