#!/bin/sh

#SBATCH --job-name="2um"
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --mail-type=ALL

module load rapids-singlecell 

# Use the specific python interpreter from rapids-singlecell
/usr/local/apps/rapids-singlecell/0.10.4/bin/python3 testRapidSiggleCell.py
