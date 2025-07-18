#!/bin/bash

# goal
# -----
# Create a liftOver chain between two reference genomes using minimap2 and associated tools.

# usage
# -----
# bash liftOverChainCreation_v2.sh $target_fasta $query_fasta $output_chain_file [$chrM_optional]
# No need of sbatch, this script is fast enough to run on an interactive job.
# [] means optional arguments. To specify an optional argument, you must specify all optional arguments to the left of
# it as well, and you must remove the square brackets for the ones you want to specify.
# target_fasta: The fasta file of the original genome assembly.
# query_fasta: The fasta file of the new genome assembly.
# output_chain_file: The output chain file that will be created.
# chrM_optional: two strings separated by a comma, the names of the mitochondrial chromosomes in the target and query
#                genomes respectively. Default is "chrM,chrM". If you specify this, we will exclude these
#                chromosomes from the chain file. WARNING: we do not checks on the user-supplied string here.

# usage of chain file
# -------------------
# load liftOver and do
# $liftOver $old_bed $output_chain_file $new_bed $unmapped
# where old_bed is a bed file in the old genome or the original genome assembly,
# output_chain_file is the chain file created by this script,
# new_bed is the output bed file in the new genome or the new genome assembly,
# and unmapped is the output file that contains the regions that could not be mapped.

# outputs
# -------
# Output file in the chain format is created at $output_chain_file.
# This is a plain text format; we are not going to describe it further here.

# stop execution if any command fails
set -e

# get current directory and go to the parent directory where the main scripts are
curr_dir=$(pwd)
cd ..

# load packages
source load_package.sh -seqkit -paf2chain -minimap2pt24

# load configuration
source config.sh

# go back to the current directory
cd "$curr_dir"

# check that the correct number of arguments were provided
if [ "$#" -lt 3 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash liftOverChainCreation_v2.sh \$target_fasta \$query_fasta \$output_chain_file [\$chrM_optional]"
    >&2 echo "For more details on what these parameters mean, see the comments in the script."
    >&2 echo "Exiting."
    exit 1;
fi

# assign arguments to variables
target_fasta=${1:-}
query_fasta=${2:-}
output_chain_file=${3:-}
chrM_optional=${4:-"chrM,chrM"}

# ensure that the fasta files exist
if [ ! -f "$target_fasta" ]; then
    >&2 echo "ERROR: target fasta file not found: $target_fasta"
    >&2 echo "Exiting."
    exit 1;
fi
if [ ! -f "$query_fasta" ]; then
    >&2 echo "ERROR: query fasta file not found: $query_fasta"
    >&2 echo "Exiting."
    exit 1;
fi

# convert paths to absolute paths
target_fasta=$(realpath "$target_fasta")
query_fasta=$(realpath "$query_fasta")
output_chain_file=$(realpath "$output_chain_file")

# ensure that the variable paf2chain is set
if [ "${paf2chain:-}" == "" ]; then
    >&2 echo "ERROR: paf2chain is not set. Please load the paf2chain package."
    >&2 echo "Exiting."
    exit 1;
fi

# set temporary directory
tmpDir="${config[scratchDir]:-}"/tmp/temp$(openssl rand -hex 6)
mkdir -p "$tmpDir"

# set temp files
tmp_target="$tmpDir"/target_filt.fasta
tmp_query="$tmpDir"/query_filt.fasta
tmp_paf="$tmpDir"/temp_paf.paf

# use seqkit to get rid of the mitochondrial chromosomes
seqkit grep -v -p "$chrM_optional" "$target_fasta" > "$tmp_target"
seqkit grep -v -p "$chrM_optional" "$query_fasta" > "$tmp_query"

# run minimap2
# this is the command recommended by the minimap2 documentation at
# https://github.com/lh3/minimap2/blob/master/cookbook.md#genome-aln
# "asm5" for sequence divergence below 1%
minimap2 -cx asm5 --cs "$tmp_target" "$tmp_query"  > "$tmp_paf"
# run paf2chain
$paf2chain -i "$tmp_paf" > "$output_chain_file"

# remove temporary directory
rm -rf "$tmpDir"