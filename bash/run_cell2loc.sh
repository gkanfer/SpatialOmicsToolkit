#!/bin/sh

#SBATCH --job-name="cell2loc"
#SBATCH --partition=gpu
#SBATCH --time=12:00:00
#SBATCH --gres=gpu:a100:4
#SBATCH --cpus-per-task=6
#SBATCH --mem=100g
#SBATCH --mail-type=ALL

module load cell2location

# Use the specific python interpreter from rapids-singlecell
python run_cell2loc.py
