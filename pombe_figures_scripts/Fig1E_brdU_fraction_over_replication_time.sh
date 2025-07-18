#!/bin/bash

# Goal: To plot BrdU fraction on individual reads versus the population median replication time.
# Usage: bash Fig1E_brdU_fraction_over_replication_time.sh /path/to/input_dir /path/to/output_dir

set -euo pipefail

# Get user-defined input/output directories
input_dir="$1"
output_dir="$2"

# Load packages
source load_package.sh -R

# Array of sample labels and corresponding BAM file prefixes
declare -A samples=(
    ["ARY017_10uM-R9"]="SacCer"
    ["ARY017_30uM-R9"]="SacCer"
    ["ARY017_100uM-R9"]="SacCer"
    ["IDS_43-R9"]="ASM294v2"
    ["IDS_47-R9"]="ASM294v2"
    ["IDS_64-R9"]="ASM294v2"
    ["IDS_49-R9"]="ASM294v2"
    ["IDS_52-R9"]="ASM294v2"
    ["IDS_65-R9"]="ASM294v2"
)

# Trep files for each species
declare -A trep_files=(
    ["SacCer"]="replication_time_s_cerevisiae.bed"
    ["ASM294v2"]="replication_time_s_pombe.bed"
)

# Get BrdU mean fraction per 1kb window per read for each sample
for sample in "${!samples[@]}"; do
    species="${samples[$sample]}"
    input_mod_bam="$input_dir/${sample}-${species}.mod.bam"
    input_trep_file="$input_dir/${trep_files[$species]}"
    output_trep_brdU_mean_per_read="$output_dir/trep_brdU_mean_per_read_${sample}-${species}.txt"

    echo "Processing $sample..."

    sbatch run_get_agg_brdU_ref_coord_int.sh "$input_mod_bam" "$input_trep_file" "$output_dir/temp_brdU_mean"
    
    grep -v '^#' "$input_trep_file" | cut -f 1-5 > "$output_dir/trep_brdU_mean_5col.txt"

    sbatch run_get_read_brdU_ref_coord_int.sh "$input_mod_bam" "$output_dir/trep_brdU_mean_5col.txt" "$output_trep_brdU_mean_per_read" "$output_dir/temp"

    rm "$output_dir/trep_brdU_mean_5col.txt"
done

# Run R script for plotting
Rscript ./plot_brdU_fraction_over_replication_time_Fig1E.R "$output_dir"
