#!/bin/sh

#SBATCH --job-name="write_xineum_fullset"
#SBATCH --cpus-per-task=20
#SBATCH --time=18:00:00
#SBATCH --mem=40g
#SBATCH --mail-type=ALL

source myconda
mamba activate squidpy-voyagerpy 

python run_xenium_fullset.py