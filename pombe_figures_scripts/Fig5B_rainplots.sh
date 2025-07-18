#!/bin/bash

# Goal: To generate rainplots (BrdU probability and fraction versus coordinates) for Fig 5B using plot_rainplot.sh
# Usage: bash Fig5B_rainplots.sh /path/to/input_dir /path/to/output_dir 

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
bam_ids_65="$input_dir/IDS_65-R9-AW2224_assembly.mod.bam"
read_id_IDS_65="97394897-18ab-4cdb-afda-87d60fe7330a"

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
out_fig5b="$output_dir/Fig5B_rainplot"


# Ensure output directories exist
mkdir -p "$out_fig5b" 

# Generate rainplots
echo "Generating rainplot for Figure 5B..."
bash "$plot_script" "$bam_ids_65" "$read_id_IDS_65" "$thymidine_window" "$out_fig5b"
echo "Figure 5B complete."

