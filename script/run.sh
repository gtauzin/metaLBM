#!/bin/bash

cd ~/metaLBM_private/
rm -rf build
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make lbm
cd ../script
./submit_job.sh
