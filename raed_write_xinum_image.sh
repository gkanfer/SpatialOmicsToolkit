#!/bin/sh

#SBATCH --job-name="write_xineum"
#SBATCH --cpus-per-task=20
#SBATCH --mem=150g
#SBATCH --mail-type=ALL

source myconda 
mamba activate squidpy-voyagerpy 

file_name="morphology_focus_0000.ome.tif"
path_morphology_focus="/data/kanferg/Sptial_Omics/playGround/Data/Xenium/Xenium_Prime_Human_Lymph_Node_Reactive_FFPE/morphology_focus"
scale=10
channels=(0 1 2 3)
prefix="morphology_focus_0000"

for ch in "${channels[@]}"
do
    pref="${prefix}_${ch}"
    python raed_write_xinum_image.py --file_name "$file_name" --path_morphology_focus "$path_morphology_focus" --scale "$scale"  --ch "$ch" --prefix "$pref"
done

file_name="morphology_focus_0001.ome.tif"
prefix="morphology_focus_0001"

for ch in "${channels[@]}"
do
    pref="${prefix}_${ch}"
    python raed_write_xinum_image.py --file_name "$file_name" --path_morphology_focus "$path_morphology_focus" --scale "$scale"  --ch "$ch" --prefix "$pref"
done

file_name="morphology_focus_0002.ome.tif"
prefix="morphology_focus_0002"

for ch in "${channels[@]}"
do
    pref="${prefix}_${ch}"
    python raed_write_xinum_image.py --file_name "$file_name" --path_morphology_focus "$path_morphology_focus" --scale "$scale"  --ch "$ch" --prefix "$pref"
done

file_name="morphology_focus_0003.ome.tif"
prefix="morphology_focus_0003"

for ch in "${channels[@]}"
do
    pref="${prefix}_${ch}"
    python raed_write_xinum_image.py --file_name "$file_name" --path_morphology_focus "$path_morphology_focus" --scale "$scale"  --ch "$ch" --prefix "$pref"
done