#!/bin/bash

# goal
# -----
# Given a plus bedgraph and a minus bedgraph and a fasta index, report (length-weighted) mean val in bedgraph in bed
# format per window of given length along the genome (also see logic flow below).
# The script name includes "genomecov" as we expect the bedgraphs to be genome-covering i.e. we want bedgraphs with
# a value at every base along the genome. But, we don't check for this, so you can violate this if you want.

# usage and inputs
#-----------------
# bash convert_genomecov_bedgraph_to_bed.sh \$plus_bedgraph \$minus_bedgraph \$fasta_fai \$interval_length \$chrM_optional
# - \$plus_bedgraph: bedgraph of plus strand data
# - \$minus_bedgraph: bedgraph of minus strand data
# - \$fasta_fai: fasta index file of genome
# - \$interval_length: length of intervals in bp to tile the genome with
# - \$chrM_optional: (default chrM) optional (regex-like) contig to exclude from analysis
#                   NOTE: this is regex-like, so chrM will match chrM, chrMito, chrMT, a_chrM etc.
#                         If you want to exclude some other contig, then you need to be careful, as
#                         for e.g. chrI will match chrI, chrII, chrIII, chrIV as all of them contain the string chrI.

# logic flow
# ----------
# - combine bedgraphs into one bed file
# - make windows along the genome in a temporary bed file
# - intersect windows with bed file
# - report (length-weighted) mean value in bedgraph in bed format per window
# NOTE: Example of length-weighted avg: let's say 3 bedgraph intervals of lengths l_1, l_2, l_3 and values v_1, v_2, v_3
#       map to a tiled genomic window, then the average reported there is
#       (l_1 * v_1 + l_2 * v_2 + l_3 * v_3) / (l_1 + l_2 + l_3).
# NOTE: Missing intervals are not counted in the average. If you want
#       to count missing intervals as 0, then you have to use bedgraphs that have 0s in the missing intervals.
# NOTE: To get a good value for the average, we do a check to ensure the mean bedgraph interval size is at least
#       10 times smaller than the tiled window size.
# NOTE: The script name includes "genomecov" as we expect the bedgraphs to be genome-covering i.e. we want bedgraphs
#       with a value at every base along the genome. But, we don't check for this, so you can violate this if you want.

# outputs
# -------
# Write in bed format to stdout.
# comments start with #, and there are no column names.
# columns are: contig, start, end, name, mean, strand, where name = "blank".

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
interval_length=$4
chrM_optional=${5:-chrM}

# check that the correct number of arguments were provided
if [ "$#" -lt 4 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash convert_genomecov_bedgraph_to_bed.sh \$plus_bedgraph \$minus_bedgraph \$fasta_fai \$interval_length \$chrM_optional"
    >&2 echo "       \$plus_bedgraph: bedgraph of plus strand data"
    >&2 echo "       \$minus_bedgraph: bedgraph of minus strand data"
    >&2 echo "       \$fasta_fai: fasta index file of genome"
    >&2 echo "       \$interval_length: length of intervals in bp to tile the genome with"
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

# check that the interval length is a positive integer
if ! [[ "$interval_length" =~ ^[0-9]+$ ]]; then
    >&2 echo "ERROR: interval length is not a positive integer"
    exit 1;
fi

# check that the interval length is greater than 0
if [ "$interval_length" -le 0 ]; then
    >&2 echo "ERROR: interval length is not greater than 0"
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

# ensure that the bed file is valid
if [ ! "$(< "$combined_bed_temp" python validate_bed_format.py --allow-float-score --never-zero-base)" == "valid" ]; then
    >&2 echo "Error: bedgraph file is not in the correct format."
    exit 1;
fi

if [ ! "$(< "$combined_bed_temp" python validate_bed_against_fai.py "$fasta_fai" )" == "valid"  ]; then
  >&2 echo "Error: bedgraph file does not have valid coordinates."
  exit 1;
fi

# ensure mean bed interval size is at least 10 times smaller than window size
combined_bed_mean_interval_size=$(awk 'BEGIN{sum=0;n=0}{sum += $3 - $2; n += 1} END{print sum/n}' "$combined_bed_temp")

if [ "$(echo "$combined_bed_mean_interval_size > $interval_length/10" | bc -l)" -eq 1 ]; then
  >&2 echo "Error: mean bed interval size is more than a tenth of the set interval length."
  exit 1;
fi

# ensure there are no self-intersections in bed file
combined_bed_temp_self_int_count=$(bedtools intersect -s -a "$combined_bed_temp" -b "$combined_bed_temp" -wao | wc -l)
combined_bed_temp_line_count=$(grep -c -E -v '^browser|^track|^#' "$combined_bed_temp")

if [ "$combined_bed_temp_self_int_count" -ne "$combined_bed_temp_line_count" ]; then
  >&2 echo "Error: bed graph file has self-intersections."
  exit 1;
fi

# make windows along the genome in a temporary bed file
windows_bed_temp="$tmpDir"/windows.bed
windows_bed_temp_no_strand="$tmpDir"/windows_no_strand.bed
bedtools makewindows -g "$fasta_fai" -w "$interval_length" > "$windows_bed_temp_no_strand"
{
  awk -v OFS="\t" '{print $1, $2, $3, "blank", 0, "+"}' "$windows_bed_temp_no_strand"
  awk -v OFS="\t" '{print $1, $2, $3, "blank", 0, "-"}' "$windows_bed_temp_no_strand"
} | grep -v "$chrM_optional" | bedtools sort -g "$fasta_fai" > "$windows_bed_temp"

# calculate mean values of bedgraph in each interval
# NOTE: we require a little over 50% of the bed graph interval to overlap with the tiled window
#       to avoid double counting i.e. to avoid the same bedgraph interval being included in two tiled windows.
# NOTE: we report a '0' instead of 'NA' or equivalent for a tiled window that has no bedgraph intervals overlapping it.
#       This is because as we have stated earlier, we expect bedgraphs to be "genome-covering", which means that there
#       must be a value at every base along the genome in the input. In spite of this, if there are missing values,
#       we assume these are zeros.
bedtools map -s -a "$windows_bed_temp" -b "$combined_bed_temp" -F 0.5000001 -c 7,8 -o sum,sum -g "$fasta_fai" -null 0 |\
  awk -v OFS="\t" '{print $1, $2, $3, $4, ($7 > 0) ? $8/$7 : 0, $6}' |\
  insert_calling_script_header "$@"

# clean up
rm -rf "$tmpDir"