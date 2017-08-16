#!/bin/bash

echo "synchronize.sh starts with argument: " $1

cd ..
filter_git=1

if [[ "$1" = "stromboli" ]]
then
   rsynct $1 ./ home/tauzin/Workspace/lbm_solver/ ${filter_git}
elif [[ "$1" = "newturb" ]]
then
    rsynct $1 ./ storage/tauzin/Workspace/lbm_solver/ ${filter_git}
elif [[ "$1" = "hydra8" ]]
then
    rsynct $1 ./ home/tauzin/Workspace/lbm_solver/ ${filter_git}
elif [[ "$1" = "galileo" ]]
then
    rsynct $1 ./ gpfs/scratch/userexternal/gtauzin0/Workspace/lbm_solver/ ${filter_git}
elif [[ "$1" = "num06" ]]
then
    rsynct $1 ./ localhdd/Workspace/lbm_solver/ ${filter_git}
fi

echo "synchronize.sh ends"
