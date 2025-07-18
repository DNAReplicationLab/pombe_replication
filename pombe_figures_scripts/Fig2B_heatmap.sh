#!/bin/bash

# Goal: To generate a heatmap with BrdU fractions for Fig 2B using plot_heatmap.sh
# Usage: bash Fig2B_heatmap.sh /path/to/input_dir /path/to/output_dir 

set -euo pipefail

# Usage message
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <input_directory> <output_directory>"
  exit 1
fi

# Set input and output directories from user arguments
input_dir="$1"
output_dir="$2"

# Set script directory
script_dir="$(dirname "$0")"
plot_script="$script_dir/plot_heatmap.sh"

# Create output subdirectory
mkdir -p "$output_dir/Fig2B_heatmap"

# Filter reads that have a fork in coordinates NC_003424.3:820000-920000
awk -F' ' '$1 == "NC_003424.3" && $2 > 820000 && $2 < 920000 {print $4}' \
  "$input_dir/left_forks-IDS_65-R9-ASM294v2.bed" \
  "$input_dir/right_forks-IDS_65-R9-ASM294v2.bed" | \
  sort | uniq > "$output_dir/Fig2B_heatmap/read_ids_Fig2B.txt"

# Obtain a bam file with the reads of interest
source load_package.sh -samtools

samtools view -h -N "$output_dir/Fig2B_heatmap/read_ids_Fig2B.txt" "$input_dir/IDS_65-R9-ASM294v2.mod.bam" > "$output_dir/Fig2B_heatmap/IDS_65-R9-ASM294v2_forks_NC_003424-3_820000-920000.bam"

samtools sort -o "$output_dir/Fig2B_heatmap/IDS_65-R9-ASM294v2_forks_NC_003424-3_820000-920000.sorted.bam" "$output_dir/Fig2B_heatmap/IDS_65-R9-ASM294v2_forks_NC_003424-3_820000-920000.bam"

samtools index "$output_dir/Fig2B_heatmap/IDS_65-R9-ASM294v2_forks_NC_003424-3_820000-920000.sorted.bam"

# Plot heatmap
bash "$plot_script" "$output_dir/Fig2B_heatmap/IDS_65-R9-ASM294v2_forks_NC_003424-3_820000-920000.sorted.bam" NC_003424.3 820000 920000 T use_ref 0.5 1000 "$output_dir/Fig2B_heatmap"

