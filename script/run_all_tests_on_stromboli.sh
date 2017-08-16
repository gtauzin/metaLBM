#!/bin/bash

echo "run_all_tests_on_stromboli starts."

echo "-- Relocating to the script folder"
cd $HOME/Workspace/lbm_solver/script/

echo "-- Synchronizing with latest stable git version"
./synchronize.sh stromboli

echo "-- SSH to stromboli to relocate to the script folder there and run lbm_stromboli.sh"
sshp stromboli "cd Workspace/lbm_solver/script/ && ./run_all_tests_stromboli.sh "

echo "run_all_tests_on_stromboli ends"
