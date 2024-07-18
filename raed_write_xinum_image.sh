#!/bin/sh

#SBATCH --job-name="write_xineum"
#SBATCH --cpus-per-task=20
#SBATCH --mem=150g
#SBATCH --mail-type=ALL

source myconda 
mamba activate squidpy-voyagerpy 
python raed_write_xinum_image.py

