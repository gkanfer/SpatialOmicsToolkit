#!/bin/sh

#SBATCH --job-name="cell2loc_reff"
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:4
#SBATCH --time=5:00:00
#SBATCH --mem=100g
#SBATCH --cpus-per-task=6
#SBATCH --mail-type=ALL

module load cell2location
python run_cell2loc_ref_genrate.py