#!/bin/bash
#
#SBATCH --job-name=lbm_2_1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --exclusive
#SBATCH --mail-user=guillaumetauzin.ut@gmail.com
#SBATCH --mail-type=ALL

export OMP_PROC_BIND=true
export OMP_PLACES=cores

cd $HOME/metaLBM_private/bin
srun ./lbm
