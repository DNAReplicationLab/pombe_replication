#!/bin/bash

# goal
# -----
# Given a trep bedgraph and a pause sensitivity bedgraph, plot the two against each other

# usage and inputs
#-----------------
# bash calculate_and_plot_trep_vs_pause_sens.sh $trep_bedgraph $pause_sens_bedgraph $fasta_fai $sequencing_summary $output_prefix
# $trep_bedgraph: bedgraph of trep values, four columns contig start end trep, space-separated, no column names,
#                 comments starting with '#'. The interval lengths must be much larger than the interval lengths
#                 in the pause sensitivity bedgraph. Otherwise, the results will be incorrect.
# $pause_sens_bedgraph: bedgraph of pause sensitivity values, four columns contig start end pause_sens, space-separated,
#                       no column names, comments starting with '#'. The interval lengths must be much smaller than
#                       the interval lengths in the trep bedgraph. Otherwise, the results will be incorrect.
# $fasta_fai: fasta index file of reference genome, produced by samtools faidx.
# $sequencing_summary: sequencing summary file produced by nanopore sequencing pipeline, to get total number of reads
#                      for normalization.
# $output_prefix: prefix for output files. if you set it to /path/to/folder/file, the output files will be
#                 /path/to/folder/file.trep_vs_pause_sens and /path/to/folder/file.trep_vs_pause_sens.png

# outputs
# -------
# We produce two files: a table and a plot.
# - The table is tab-separated with five columns with column names and comments
#   starting with '#'. The columns are: contig, start, end, trep, sum_sens.
#   First four columns are a repeat of the input trep bedgraph. The fifth column
#   is the sum of pause sensitivity values of all intervals overlapping with the trep interval.
# - The plot is a scatter plot of trep vs sum_sens. The x-axis is trep and the y-axis is sum_sens.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -bedtools -R -python -miller

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
t_rep=${1:-}
pause_sens=${2:-}
fasta_fai=${3:-}
sequencing_summary=${4:-}
output_prefix=${5:-}

# check that the correct number of arguments were provided
if [ "$#" -ne 5 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash calculate_and_plot_trep_vs_pause_sens.sh \$trep_bedgraph \$pause_sens_bedgraph \$fasta_fai \$sequencing_summary \$output_prefix"
    >&2 echo "       \$trep_bedgraph: bedgraph of trep values, four columns contig start end trep"
    >&2 echo "       \$pause_sens_bedgraph: bedgraph of pause sensitivity values, four columns contig start end pause_sens"
    >&2 echo "       \$fasta_fai: fasta index file of reference genome, produced by samtools faidx"
    >&2 echo "       \$sequencing_summary: sequencing summary file produced by nanopore sequencing pipeline"
    >&2 echo "       \$output_prefix: prefix for output files."
    >&2 echo "For more information, see the header of this script."
    exit 1;
fi

# check that the input files exist
if [ ! -f "$t_rep" ]; then
    >&2 echo "ERROR: trep bedgraph file $t_rep does not exist"
    exit 1;
fi
if [ ! -f "$pause_sens" ]; then
    >&2 echo "ERROR: pause sensitivity bedgraph file $pause_sens does not exist"
    exit 1;
fi
if [ ! -f "$fasta_fai" ]; then
    >&2 echo "ERROR: fasta index file $fasta_fai does not exist"
    exit 1;
fi
if [ ! -f "$sequencing_summary" ]; then
    >&2 echo "ERROR: sequencing summary file $sequencing_summary does not exist"
    exit 1;
fi

# check that output prefix is not blank
if [ -z "$output_prefix" ]; then
    >&2 echo "ERROR: output prefix is blank"
    exit 1;
fi

# make the output directory if it does not exist
output_dir=$(dirname "$output_prefix")
mkdir -p "$output_dir"

# make temporary files
t_rep_temp=$(mktemp -p "$tmpDir" trep.XXXXXX)
pause_sens_temp=$(mktemp -p "$tmpDir" pause_sens.XXXXXX)

# get total number of reads first
n_reads=$(mlr --tsv --skip-comments --headerless-csv-output count "$sequencing_summary")

# convert bedgraphs to bed files
awk -v OFS="\t" '!/^#/ {if($2 != $3){print $1, $2, $3, $4}}' "$t_rep" | bedtools sort -g "$fasta_fai" > "$t_rep_temp"
awk -v OFS="\t" '!/^#/ {if($2 != $3){print $1, $2, $3, $4}}' "$pause_sens" | bedtools sort -g "$fasta_fai" \
  > "$pause_sens_temp"

# ensure that the bed files are valid
if [ ! "$(< "$t_rep_temp" python validate_bed_format.py --allow-float-score)" == "valid" ]; then
    >&2 echo "Error: trep file is not in the correct format."
    exit 1;
fi

if [ ! "$(< "$pause_sens_temp" python validate_bed_format.py --allow-float-score)" == "valid" ]; then
    >&2 echo "Error: pause sensitivity file is not in the correct format."
    exit 1;
fi

if [ ! "$(< "$t_rep_temp" python validate_bed_against_fai.py "$fasta_fai" )" == "valid"  ]; then
  >&2 echo "Error: trep file does not have valid coordinates."
  exit 1;
fi

if [ ! "$(< "$pause_sens_temp" python validate_bed_against_fai.py "$fasta_fai" )" == "valid"  ]; then
  >&2 echo "Error: pause sensitivity file does not have valid coordinates."
  exit 1;
fi

# ensure mean pause sensitivity interval size is at least 10 times smaller than trep interval size
pause_sens_mean_interval_size=$(awk 'BEGIN{sum=0;n=0}{sum += $3 - $2; n += 1} END{print sum/n}' "$pause_sens_temp")
t_rep_mean_interval_size=$(awk 'BEGIN{sum=0;n=0}{sum += $3 - $2; n += 1} END{print sum/n}' "$t_rep_temp")

if [ "$(echo "$pause_sens_mean_interval_size > $t_rep_mean_interval_size/10" | bc -l)" -eq 1 ]; then
  >&2 echo "Error: mean pause sensitivity interval size is more than a tenth of the mean trep interval size."
  exit 1;
fi

# ensure there are no self-intersections in either bed file
pause_sens_self_int_count=$(bedtools intersect -a "$pause_sens_temp" -b "$pause_sens_temp" -wao | wc -l)
t_rep_self_int_count=$(bedtools intersect -a "$t_rep_temp" -b "$t_rep_temp" -wao | wc -l)
pause_sens_line_count=$(grep -c -E -v '^browser|^track|^#' "$pause_sens_temp")
t_rep_line_count=$(grep -c -E -v '^browser|^track|^#' "$t_rep_temp")

if [ "$pause_sens_self_int_count" -ne "$pause_sens_line_count" ]; then
  >&2 echo "Error: pause sensitivity file has self-intersections."
  exit 1;
fi

if [ "$t_rep_self_int_count" -ne "$t_rep_line_count" ]; then
  >&2 echo "Error: trep file has self-intersections."
  exit 1;
fi

# calculate the sum of pause sensitivity values for each trep interval
# and normalize summed sensitivity to per million reads
bedtools map -a "$t_rep_temp" -b "$pause_sens_temp" -c 4 -o sum -g "$fasta_fai" -null 0 |\
  awk -v n_reads="$n_reads" -v OFS="\t" '{print $1, $2, $3, $4, $5/n_reads*1000000}' |\
  sed '1icontig\tstart\tend\ttrep\tsum_sens' |\
  sed '1i# sensitivity normalized to per million reads, total num reads: '"$n_reads"'' |\
  insert_calling_script_header "$@" >  "$output_prefix".trep_vs_pause_sens

# make the scatter plot
cd plotting_and_short_analyses;
< "$output_prefix".trep_vs_pause_sens Rscript plot_scatter.R trep sum_sens\
  "$output_prefix".trep_vs_pause_sens.png "Median Replication time (min)" "Pause sensitivity per million reads"\
  15,65 0,auto 0.1 1 0 5

# clean up
rm -rf "$tmpDir"