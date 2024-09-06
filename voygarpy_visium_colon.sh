#!/bin/sh

#SBATCH --job-name="vp_full_analysis"
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:2
# #SBATCH --time=5:00:00
#SBATCH --mem=50g
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=ALL

source myconda
mamba activate squidpy-voyagerpy

/data/kanferg/conda/envs/squidpy-voyagerpy/bin/python3 voygarpy_visium_colon.py