#!/bin/sh

#SBATCH --job-name="stlearn LR analysis"
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --time=5:00:00
#SBATCH --mem=50g
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL

source myconda
mamba activate stlearn-env

/gpfs/gsfs10/users/kanferg/conda/envs/stlearn-env/bin/python run_stlearn_LR_analysis_visium.py
