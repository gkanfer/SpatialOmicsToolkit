#!/bin/sh

#SBATCH --job-name="write_xineum_smallset"
#SBATCH --cpus-per-task=20
#SBATCH --time=8:00:00
#SBATCH --mem=280g
#SBATCH --mail-type=ALL

source myconda
mamba activate squidpy-voyagerpy 

python run_xenium_smallSet.py
