#!/bin/bash

# goal
# -----
# Given a plus bedgraph and a minus bedgraph, report (length-weighted) mean val per bed entry
# (also see logic flow below).

# usage and inputs
#-----------------
# bash calculate_bed_overlap_bedgraph_signal_per_bed_entry.sh \$plus_bedgraph \$minus_bedgraph \$fasta_fai \$bed_file \$chrM_optional
# - \$plus_bedgraph: bedgraph of plus strand data
# - \$minus_bedgraph: bedgraph of minus strand data
# - \$fasta_fai: fasta index file of genome
# - \$bed_file: bed file to intersect with
# - \$chrM_optional: (default chrM) optional (regex-like) contig to exclude from analysis
#                   NOTE: this is regex-like, so chrM will match chrM, chrMito, chrMT, a_chrM etc.
#                         If you want to exclude some other contig, then you need to be careful, as
#                         for e.g. chrI will match chrI, chrII, chrIII, chrIV as all of them contain the string chrI.

# logic flow
# ----------
# - combine bedgraphs into one bed file
# - intersect combined bed file with given bed file in a stranded manner
# - report (length-weighted) mean value in bedgraph per bed entry
# NOTE: Example of length-weighted avg: let's say 3 bedgraph intervals of lengths l_1, l_2, l_3 and values v_1, v_2, v_3
#       map to a bed entry, then the average reported there is
#       (l_1 * v_1 + l_2 * v_2 + l_3 * v_3) / (l_1 + l_2 + l_3).
# NOTE: Missing intervals are not counted in the average. If you want
#       to count missing intervals as 0, then you have to use bedgraphs that have 0s in the missing intervals.
# NOTE: While it is desirable that each bed entry covers several bedgraph intervals so that our average
#       is not too noisy, we do not explicitly do such a check here.

# outputs
# -------
# Write in bed format to stdout, adding an additional column to the input bed file which is the (length-weighted) mean.
# comments start with #, and there are no column names.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages
source load_package.sh -bedtools -python

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
plus_bedgraph=$1
minus_bedgraph=$2
fasta_fai=$3
bed_file=$4
chrM_optional=${5:-chrM}

# check that the correct number of arguments were provided
if [ "$#" -lt 4 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash $0 \$plus_bedgraph \$minus_bedgraph \$fasta_fai \$bed_file \$chrM_optional"
    >&2 echo "       \$plus_bedgraph: bedgraph of plus strand data"
    >&2 echo "       \$minus_bedgraph: bedgraph of minus strand data"
    >&2 echo "       \$fasta_fai: fasta index file of genome"
    >&2 echo "       \$bed_file: bed file to intersect with"
    >&2 echo "       \$chrM_optional: (default chrM) optional (regex-like) contig to exclude from analysis"
    >&2 echo "                        NOTE: this is regex-like, so for e.g. chrI will match chrI, chrII, chrIII, etc."
    >&2 echo "NOTE: It is recommended that you read the header of this script to understand what it does."
    exit 1;
fi

# check that the input files exist
if [ ! -f "$plus_bedgraph" ]; then
    >&2 echo "ERROR: plus bedgraph file does not exist"
    exit 1;
fi

if [ ! -f "$minus_bedgraph" ]; then
    >&2 echo "ERROR: minus bedgraph file does not exist"
    exit 1;
fi

if [ ! -f "$fasta_fai" ]; then
    >&2 echo "ERROR: fasta index file does not exist"
    exit 1;
fi

if [ ! -f "$bed_file" ]; then
    >&2 echo "ERROR: bed file does not exist"
    exit 1;
fi

# make a temporary bed file
combined_bed_temp="$tmpDir"/combined.bed

# convert bedgraphs to bed files. we add two more columns: length of interval, and signal * length of interval.
# these will be used to calculate the weighted mean signal per interval.
{
  awk -v OFS="\t" '!/^#/ {print $1, $2, $3, "blank", $4, "+", $3 - $2, ($3 - $2) * $4}' "$plus_bedgraph"
  awk -v OFS="\t" '!/^#/ {print $1, $2, $3, "blank", $4, "-", $3 - $2, ($3 - $2) * $4}' "$minus_bedgraph"
} | grep -v "$chrM_optional" | bedtools sort -g "$fasta_fai" > "$combined_bed_temp"

# ensure that the combined bed file is valid
if [ ! "$(< "$combined_bed_temp" python validate_bed_format.py --allow-float-score --never-zero-base)" == "valid" ]; then
    >&2 echo "Error: bedgraph file is not in the correct format."
    exit 1;
fi

if [ ! "$(< "$combined_bed_temp" python validate_bed_against_fai.py "$fasta_fai" )" == "valid"  ]; then
  >&2 echo "Error: bedgraph file does not have valid coordinates."
  exit 1;
fi

# ensure that the input bed file is valid
if [ ! "$(< "$bed_file" python validate_bed_format.py --allow-float-score --six-columns)" == "valid" ]; then
    >&2 echo "Error: bed file is not in the correct format."
    exit 1;
fi

if [ ! "$(< "$bed_file" python validate_bed_against_fai.py "$fasta_fai" )" == "valid"  ]; then
  >&2 echo "Error: bed file does not have valid coordinates."
  exit 1;
fi

# ensure there are no self-intersections in the combined bed file
combined_bed_temp_self_int_count=$(bedtools intersect -s -a "$combined_bed_temp" -b "$combined_bed_temp" -wao | wc -l)
combined_bed_temp_line_count=$(grep -c -E -v '^browser|^track|^#' "$combined_bed_temp")

if [ "$combined_bed_temp_self_int_count" -ne "$combined_bed_temp_line_count" ]; then
  >&2 echo "Error: bed graph file has self-intersections."
  exit 1;
fi

# make a temporary bed file and sort the input file
bed_file_temp="$tmpDir"/bed_file.bed
grep -v "$chrM_optional" "$bed_file" | bedtools sort -g "$fasta_fai" > "$bed_file_temp"

# calculate mean values of bedgraph in each interval
# NOTE: we require a little over 50% of the bed graph interval to overlap with the bed window.
# use awk to print all columns except last two and then add a new column that is the ratio of the last two columns
bedtools map -s -a "$bed_file_temp" -b "$combined_bed_temp" -F 0.5000001 -c 7,8 -o sum,sum -g "$fasta_fai" -null 0 |\
  awk '{for (i=1;i<NF-1;i++) printf "%s\t",$i; printf "%s\n", ($((NF-1)) > 0) ? $NF/$((NF-1)) : "NA"}' |\
  insert_calling_script_header "$@"

# clean up
rm -rf "$tmpDir"