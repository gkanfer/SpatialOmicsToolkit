#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --gres=gpu:p100:1
#SBATCH --job-name=scvitools_test

module load scvitoos/1.0.4.gpu
set -e
cd /data/$USER/scvitools_test/
python-scvitools run_scvi.pys