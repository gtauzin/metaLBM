#!/bin/bash
#
#SBATCH --job-name=lbm_1_1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=gpu-out.%j
#SBATCH --error=gpu-err.%j
#SBATCH --time=00:05:00
#SBATCH --partition=gpus
#SBATCH --mail-user=guillaumetauzin.ut@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH --gres=gpu:1

export OMP_PROC_BIND=true
export OMP_PLACES=cores

cd $HOME/metaLBM_private/bin
srun ./lbm
