#!/bin/bash

# Goal: Detect replication pauses for Fig 3C using run_rDNA_detectSummary.sh
# Usage: bash Fig3C_pause_detection.sh /path/to/input_dir /path/to/output_dir 

set -euo pipefail

# Check input arguments
if [[ $# -ne 2 ]]; then
  echo "Usage: $0 /path/to/input_dir /path/to/output_dir"
  exit 1
fi

# Define input and output directories
input_dir="$1"
output_dir="$2"

# Define file paths and parameters
mod_bam="$input_dir/IDS_65-R9-rDNA.mod.bam"
output_directory="$output_dir/Fig3C"
ref_fasta="$input_dir/reference_rDNA.fasta"
bed_file_with_features="$input_dir/rDNA_reference_regions_on_sp_rdna_1_repeat.bed"
pos_boundary=6772
number=25

# Create output directory
mkdir -p "$output_directory"

# Submit SLURM job with all arguments in a single command line
sbatch run_rDNA_detectSummary.sh "$mod_bam" "$output_directory" "$ref_fasta" "" "" "$bed_file_with_features" "" "$number" "$pos_boundary" ""
