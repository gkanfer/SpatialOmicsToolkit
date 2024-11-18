#!/bin/sh

#SBATCH --job-name="scvi_env analysis"
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:3
#SBATCH --time=24:00:00
#SBATCH --mem=6g
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL

source myconda
mamba activate scvi_env

/gpfs/gsfs10/users/kanferg/conda/envs/scvi_env/bin/python Breast_Cancer_scVI_LabelTransfer.py

