#!/bin/bash

# Goal: To generate rainplots (BrdU probability and fraction versus coordinates) for Fig 1B–1D using plot_rainplot.sh
# Usage: bash Fig1B-D_rainplots.sh /path/to/input_dir /path/to/output_dir 

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
bam_ids_47="$input_dir/IDS_47-R9-ASM294v2.mod.bam"
read_id_47="1dc9d397-aeea-4d1a-8b02-341d5056f4bc"

bam_ary017="$input_dir/ARY017_30uM-R9-SacCer.mod.bam"
read_id_ary017="adaa3d74-af88-47ba-96ff-5781310d9fba"

bam_ids_65="$input_dir/IDS_65-R9-ASM294v2.mod.bam"
read_id_ids_65="f0c17403-b7ab-40a2-b98a-f2df2d868ebf"

# Validate input files
for bam_file in "$bam_ids_47" "$bam_ary017" "$bam_ids_65"; do
  if [ ! -f "$bam_file" ]; then
    echo "Error: Input BAM file not found: $bam_file"
    exit 1
  fi
done

# Parameters
thymidine_window="300"

# Output directories
out_fig1b="$output_dir/Fig1B_rainplot"
out_fig1c="$output_dir/Fig1C_rainplot"
out_fig1d="$output_dir/Fig1D_rainplot"

# Ensure output directories exist
mkdir -p "$out_fig1b" "$out_fig1c" "$out_fig1d"

# Generate rainplots
echo "Generating rainplot for Figure 1B..."
bash "$plot_script" "$bam_ids_47" "$read_id_47" "$thymidine_window" "$out_fig1b"
echo "Figure 1B complete."

echo "Generating rainplot for Figure 1C..."
bash "$plot_script" "$bam_ary017" "$read_id_ary017" "$thymidine_window" "$out_fig1c"
echo "Figure 1C complete."

echo "Generating rainplot for Figure 1D..."
bash "$plot_script" "$bam_ids_65" "$read_id_ids_65" "$thymidine_window" "$out_fig1d"
echo "Figure 1D complete."
