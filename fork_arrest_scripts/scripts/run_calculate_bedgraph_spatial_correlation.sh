#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J bedgSpatCorr
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# goal
# -----
# Calculate spatial correlation function from bedgraphs

# usage
#------
# sbatch run_calculate_bedgraph_spatial_correlation.sh bedgraph_file max_distance fasta_fai
# Can use bash but running time might be a couple of hours or so depending on the input files/parameters
# bedgraph_file: bedgraph file containing pause sensitivities.
# max_distance: maximum distance to consider for nearest neighbour distance statistics
# fasta_fai: fasta index file for the genome, usually ends in .fasta.fai. or .fa.fai.
# NOTE: for more details about the first two parameters, see the script invoked by this script

# outputs
# -------
# A tab-separated plain-text table to stdout. For more details, see the script invoked by this script.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -python

# load configuration
source config.sh

# function to set calling script information
insert_calling_script_header() {
  sed '1i'\
'# from commit '"${COMMITSTR:-NA}"' generated at '"${TIMENOW:-NA}"' by '"${config[name]:-NA}"' <'"${config[email]:-NA}"'>\n'\
'# script: '"$0"'\n'\
'# arguments: '"$*"'\n'\
"# slurm job name: ${SLURM_JOB_NAME:-NA}"
}

# assign arguments to variables
bedgraph_file="$1"
max_distance="$2"
fasta_fai="$3"

# check that the correct number of arguments were provided
if [ "$#" -lt 3 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: sbatch run_calculate_theoretical_pause_neighbour_distance_from_sensitivity.sh bedgraph_file max_distance fasta_fai"
    >&2 echo "For more details on what these parameters mean, see the comments in the script."
    exit 1;
fi

# check that the input files exist
if [ ! -f "$bedgraph_file" ]; then
    >&2 echo "ERROR: input file $bedgraph_file does not exist"
    exit 1;
fi

if [ ! -f "$fasta_fai" ]; then
    >&2 echo "ERROR: input file $fasta_fai does not exist"
    exit 1;
fi

# check that max_distance is a positive integer
if ! [[ "$max_distance" =~ ^[0-9]+$ ]]; then
    >&2 echo "ERROR: max_distance must be a positive integer"
    exit 1;
fi

# check that max_distance is greater than 0
if [ "$max_distance" -lt 1 ]; then
    >&2 echo "ERROR: max_distance must be greater than 0"
    exit 1;
fi

# validate the sensitivity bedgraph
# NOTE: we use a big coverage tolerance below because some bedgraphs have signals over regions
# like the tRNAs, rDNA regions removed. We've set it to roughly the total size of these regions in S. cerevisiae.
validity=$(< "$bedgraph_file" grep -v -E '^browser|^track|^#' | awk 'BEGIN{OFS="\t"}{print $1, $2, $3}' |\
  python validate_bed_against_fai.py --check-genome-cov-to-tol 30000  --check-no-overlap "$fasta_fai")

if [ "$validity" != "valid" ]; then
    >&2 echo "ERROR: sensitivity file $bedgraph_file is not valid"
    exit 1;
fi

# run the script to calculate pause distance statistics and send the output to stdout
< "$bedgraph_file" python calculate_bedgraph_spatial_autocorrelation.py "$max_distance" |\
 insert_calling_script_header "$@"