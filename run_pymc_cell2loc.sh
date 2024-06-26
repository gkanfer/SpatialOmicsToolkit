#!/bin/sh

#SBATCH --job-name="pymc_cell2loc"
#SBATCH --partition=gpu
#SBATCH --time=12:00:00
#SBATCH --gres=gpu:a100:4
#SBATCH --cpus-per-task=12
#SBATCH --mem=120g
#SBATCH --mail-type=ALL

module load pymc/5
python-pymc run_pymc_cell2loc.py
