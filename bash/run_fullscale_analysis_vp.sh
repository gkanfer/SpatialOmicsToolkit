#!/bin/sh

#SBATCH --job-name="vp_full_analysis"
#SBATCH --time=6:00:00
#SBATCH --mem=150g
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=ALL

source myconda
mamba activate squidpy-voyagerpy

/data/kanferg/conda/envs/squidpy-voyagerpy/bin/python3 run_fullscale_analysis_vp.py