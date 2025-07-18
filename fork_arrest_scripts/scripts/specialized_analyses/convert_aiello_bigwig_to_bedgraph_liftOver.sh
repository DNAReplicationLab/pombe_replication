#!/bin/bash
#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH -p ei-short
#SBATCH -J convertBwLiftOver
#SBATCH --mail-type=END,FAIL
#SBATCH --time=00:59:59
#SBATCH --constraint=""

# goal
# ----
# Convert a bigwig file to bedgraph format and lift over the coordinates to a new genome, merging values that map to the same region.


# usage
# -----
# bash convert_bigwig_liftover.sh <input_bigwig> <liftover_chain> <new_genome_fai> <output_bedgraph>
# input_bigwig: bigwig file to convert
# liftover_chain: chain file to use for liftover
# new_genome_fai: fai file for the new genome
# output_bedgraph: output bedgraph file

# outputs
# -------
# output_bedgraph: bedgraph file with the converted values

# stop execution if any step fails
set -e

# load bigwigToBedGraph and liftOver
cd ..;
source load_package.sh -bigWigToBedGraph -liftOver

# load git repo labels
source load_git_repo_labels.sh

# load configuration
source config.sh

# set temporary directory
tmpDir="${config[scratchDir]:-}"/tmp/temp$(openssl rand -hex 6)
mkdir -p "$tmpDir"

# Input arguments
input_bigwig=${1:-}
liftover_chain=${2:-}
new_genome_fai=${3:-}
output_bedgraph=${4:-}

# check that the correct number of arguments were provided
if [ "$#" -lt 4 ]; then
  >&2 echo "Error: 4 arguments are required"
  >&2 echo "bash convert_bigwig_liftover.sh <input_bigwig> <liftover_chain> <new_genome_fai> <output_bedgraph>"
  >&2 echo "input_bigwig: bigwig file to convert"
  >&2 echo "liftover_chain: chain file to use for liftover" 
  >&2 echo "new_genome_fai: fai file for the new genome"
  >&2 echo "output_bedgraph: output bedgraph file"
  exit 1
fi

# check that the three input files exist
if [ ! -f "$input_bigwig" ]; then
  >&2 echo "Error: input_bigwig file not found: $input_bigwig"
  exit 1
fi

if [ ! -f "$liftover_chain" ]; then
  >&2 echo "Error: liftover_chain file not found: $liftover_chain"
  exit 1
fi

if [ ! -f "$new_genome_fai" ]; then
  >&2 echo "Error: new_genome_fai file not found: $new_genome_fai"
  exit 1
fi

# make output folder if needed
mkdir -p "$(dirname "$output_bedgraph")"

# function to set calling script information
insert_calling_script_header() {
  sed '1i'\
'# from commit '"${COMMITSTR:-NA}"' generated at '"${TIMENOW:-NA}"' by '"${config[name]:-NA}"' <'"${config[email]:-NA}"'>\n'\
'# script: '"$0"'\n'\
'# arguments: '"$*"'\n'\
"# slurm job name: ${SLURM_JOB_NAME:-NA}"
}

# Temporary files
temp_bedgraph=$(mktemp --tmpdir="$tmpDir")
temp_bed=$(mktemp --tmpdir="$tmpDir")
temp_lifted_bed=$(mktemp --tmpdir="$tmpDir")
temp_unmapped_bed=$(mktemp --tmpdir="$tmpDir")

# Convert bigwig to bedgraph
echo "Converting bigwig to bedgraph..."
bigWigToBedGraph "$input_bigwig" "$temp_bedgraph"

# Convert bedgraph to bed
echo "Converting bedgraph to bed..."
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "name", $4, "+"}' "$temp_bedgraph" > "$temp_bed"

# Lift over the bed file
echo "Lifting over the bed file..."
$liftOver "$temp_bed" "$liftover_chain" "$temp_lifted_bed" "$temp_unmapped_bed"

# convert the lifted bed back to a bedgraph and merge overlapping intervals
echo "Converting lifted bed to bedgraph..."
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $5}' "$temp_lifted_bed"  |\
    bash merge_overlapping_bedgraph_intervals.sh "$new_genome_fai" |\
        insert_calling_script_header "$@" >  "$output_bedgraph"

# Clean up temporary files
rm -rf "$tmpDir"

# Print completion message
echo "Conversion complete. Output bedgraph: $output_bedgraph"
echo "NOTE: due to loss from liftOver etc., only use these files for visualization"
echo "There are other scripts here which do a better job in calculations"