#!/bin/bash

# Goal: To generate rainplots (BrdU probability and fraction versus coordinates) for Fig 2A using plot_rainplot.sh
# Usage: bash Fig2A_rainplot.sh /path/to/input_dir /path/to/output_dir 

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

# Input BAM file and read ID
bam_ids_65="$input_dir/IDS_65-R9-ASM294v2.mod.bam"
read_id_ids_65="30c378af-16c8-4cb8-9bc5-c4421350c0bc"

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
out_fig2a="$output_dir/Fig2A_rainplot"

# Ensure output directories exist
mkdir -p "$out_fig2a"

# Generate rainplot
echo "Generating rainplot for Figure 2A..."
bash "$plot_script" "$bam_ids_65" "$read_id_ids_65" "$thymidine_window" "$out_fig2a"
echo "Figure 2A complete."
