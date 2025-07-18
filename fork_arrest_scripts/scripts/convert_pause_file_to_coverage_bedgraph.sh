#!/bin/bash

# goal
# -----
# Given a pause file, extract valid pauses, grow by +- some bp, tile genome, and send pauses per interval to bedgraph.

# usage
#------
# bash convert_pause_file_to_coverage_bedgraph.sh pause_file slop_size genome_window_size fasta_fai
# pause_file: pause file in our usual tab-separated format. Not gonna go into details of the format here.
#             See convert_pause_file_to_bed.py or other scripts for details on pause format.
# slop_size: number of bp to grow pauses by. If 0, then no growing, and output 1bp-sized intervals.
#            use this to express uncertainty in pause location.
# genome_window_size: size of window to tile genome with. Choose to be much larger than slop_size.
# fasta_fai: fasta index file, ends in .fai. Used to get contig lengths.

# outputs
# -------
# Send to stdout a bedgraph file with the following columns (space-separated, no column name, comments start with #):
# contig, start, end, pause_coverage

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -python -miller

# load configuration
source config.sh

# set temporary directory
tmpDir="${config[scratchDir]:-}"/tmp/temp$(openssl rand -hex 6)
mkdir -p "$tmpDir"

# function to set calling script information
insert_calling_script_header() {
  sed '1i'\
'# from commit '"${COMMITSTR:-NA}"' generated at '"${TIMENOW:-NA}"' by '"${config[name]:-NA}"' <'"${config[email]:-NA}"'>\n'\
'# script: '"$0"'\n'\
'# arguments: '"$*"'\n'\
"# slurm job name: ${SLURM_JOB_NAME:-NA}"
}

# assign arguments to variables
pause_file=${1:-}
slop_size=${2:-}
genome_window_size=${3:-}
ref=${4:-}

# check that the correct number of arguments were provided
if [ "$#" -lt 1 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash convert_pause_file_to_coverage_bedgraph.sh pause_file slop_size genome_window_size ref"
    >&2 echo "pause_file: pause file in our usual tab-separated format. Not gonna go into details of the format here."
    >&2 echo "slop_size: number of bp to grow pauses by. If 0, then no growing, and output 1bp-sized intervals."
    >&2 echo "genome_window_size: size of window in bp to tile genome with. Choose to be much smaller than slop_size."
    >&2 echo "ref: fasta index file, ends in .fai. Used to get contig lengths."
    exit 1;
fi

# check that the input files exist
if [ ! -f "$pause_file" ]; then
    >&2 echo "ERROR: pause_file $pause_file does not exist"
    exit 1;
fi

if [ ! -f "$ref" ]; then
    >&2 echo "ERROR: ref $ref does not exist"
    exit 1;
fi

# check that slop size is larger than genome window size
if [ "$(echo "$slop_size >= 5 * $genome_window_size" | bc -l)" -eq 0 ]; then
    >&2 echo "ERROR: slop_size $slop_size needs to be at least 5 times genome_window_size $genome_window_size"
    exit 1;
fi

# make temporary files
tmp_pause_bed_file=$(mktemp -p "$tmpDir" tmp_pause_bed_file.XXXXXX)
tmp_bed_file=$(mktemp -p "$tmpDir" tmp_bed_file.XXXXXX)

# convert pause file to bed
< "$pause_file" python convert_pause_file_to_bed.py --discardUnkeptPauses  > "$tmp_pause_bed_file"

# count number of pauses
pause_count=$(mlr --tsv --skip-comments --headerless-csv-output count "$tmp_pause_bed_file")

# grow pauses by slop_size
bedtools slop -i "$tmp_pause_bed_file" -g "$ref" -b "$slop_size" > "$tmp_bed_file"

# get coverage per interval and normalize to pause count and output to bedgraph
bash convert_bed_to_coverage.sh "$ref" "$genome_window_size" "$tmp_bed_file" |\
  bash normalize_bedgraph_to_sum.sh "$pause_count" |\
  insert_calling_script_header "$@"

# remove temporary files
rm -rf "$tmpDir"