#!/bin/sh

#SBATCH --job-name="cell2loc_build"
# #SBATCH --partition=gpu
# #SBATCH --gres=gpu:a100:2
#SBATCH --cpus-per-task=12
#SBATCH --mem=120g
#SBATCH --mail-type=ALL

module load cell2location
python pymc_data_gnr.py


