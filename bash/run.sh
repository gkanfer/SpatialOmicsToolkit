#!/bin/sh

#SBATCH --job-name="8um"
# SBATCH --partition=gpu
# #SBATCH --gres=gpu:a100:1
# #SBATCH --gres=gpu:a100:1
# #SBATCH --time=3:00:00
#SBATCH --mem=100g
#SBATCH --cpus-per-task=6
#SBATCH --mail-type=ALL

# module load rapids-singlecell 

# python test.py
source myconda
mamba activate squidpy

python test_008.py
