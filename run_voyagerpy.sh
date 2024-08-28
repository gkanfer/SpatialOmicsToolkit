#!/bin/sh

#SBATCH --job-name="voygerpy"
# #SBATCH --exclusive
#SBATCH --time=6:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:v100x:1
#SBATCH --mem=373g
#SBATCH --cpus-per-task=36
#SBATCH --mail-type=ALL

source myconda 
mamba activate squidpy-voyagerpy

/data/kanferg/conda/envs/squidpy-voyagerpy/bin/python3 run_voyagerpy.py
# /data/kanferg/conda/envs/squidpy-voyagerpy/bin/python3 run_voyagerpy_colon.py