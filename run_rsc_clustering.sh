#!/bin/sh

#SBATCH --job-name="rapids-singlecell"
#SBATCH --partition=gpu
#SBATCH --time=3:00:00
#SBATCH --gres=gpu:a100:4
#SBATCH --cpus-per-task=10
#SBATCH --mem=50g
#SBATCH --mail-type=ALL

module load rapids-singlecell
python3-rsc run_rsc_clustering.py