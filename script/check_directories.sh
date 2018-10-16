#!/bin/bash

cd ../$1
mkdir -p ../output ../output/outputVTR ../output/outputBackup ../output/outputPlot ../output/outputAnimation> /dev/null
mkdir -p ../input/inputJSON > /dev/null
mkdir -p build bin log docs > /dev/null
