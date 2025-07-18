#!/bin/bash

# Goal: To generate rainplots (BrdU probability and fraction versus coordinates) for Fig 3B using plot_rainplot.sh
# Usage: bash Fig3B_rainplots.sh /path/to/input_dir /path/to/output_dir 

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
plot_script="$script_dir/plot_rainplot.sh"

# Input BAM files and read IDs
bam_ids_65="$input_dir/IDS_65-R9-ASM294v2.mod.bam"
read_id_IDS_65_top="84181241-826d-47c2-b75a-71b32ba820b4"
read_id_IDS_65_bottom="782229f5-14c5-4194-8a8d-9f73cdaa4c69"

# Validate input files
for bam_file in "$bam_ids_65"; do
  if [ ! -f "$bam_file" ]; then
    echo "Error: Input BAM file not found: $bam_file"
    exit 1
  fi
done

# Parameters
thymidine_window="300"

# Output directories
out_fig3b_top="$output_dir/Fig3B_top_rainplot"
out_fig3b_bottom="$output_dir/Fig3B_bottom_rainplot"

# Ensure output directories exist
mkdir -p "$out_fig3b_top" "$out_fig3b_bottom" 

# Generate rainplots
echo "Generating rainplot for Figure 3B top..."
bash "$plot_script" "$bam_ids_65" "$read_id_IDS_65_top" "$thymidine_window" "$out_fig3b_top"
echo "Figure 3B top complete."

echo "Generating rainplot for Figure 3B bottom..."
bash "$plot_script" "$bam_ids_65" "$read_id_IDS_65_bottom" "$thymidine_window" "$out_fig3b_bottom"
echo "Figure 3B bottom complete."

