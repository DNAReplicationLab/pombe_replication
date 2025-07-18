#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 10
#SBATCH -p ei-medium
#SBATCH -J calcBrduAmtVsGenomeWindowOnFeatures
#SBATCH --mail-type=END,FAIL
#SBATCH --time 23:59:59
#SBATCH --constraint=""

# This script relies on calculate_brdu_amount_vs_genome_window_on_forks.sh.
# Basically, we act like the feature file provided here is a right fork file.
# Please read the comments in that script for more information.

bash calculate_brdu_amount_vs_genome_window_on_forks.sh "$1" "$2" "$3" /dev/null "$4";