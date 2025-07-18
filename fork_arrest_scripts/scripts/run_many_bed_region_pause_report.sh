#!/bin/bash

# goal
# -----
# Given a dataset with detected pauses and a list of bed files, get a pause report for each bed file,
# once for all forks and then splitting by head-on and co-directional forks.

# usage
#------
# bash run_many_bed_region_pause_report.sh input_file.json
# - no need of sbatch as we only perform simple steps before launching many slurm jobs here.
# - input_file.json is a json file with the following format:
#   NOTE: fields can be in any order and any additional fields are ignored.
#   NOTE: be careful to not have any typos in the json file.
# {
#   "mod_bam": "/path/to/mod_bam_file.bam",
#   "forksense_dir": "/path/to/fork_sense_directory",
#   "pause_file": "/path/to/pause_file",
#   "op_dir": "/path/to/output_directory",
#   "dataset": "dataset1",
#   "analysis_label": "12jun23",
#   "mod_bam_left": "/path/to/mod_bam_left_file.bam",
#   "mod_bam_right": "/path/to/mod_bam_right_file.bam",
#   "fasta": "/path/to/fasta_file.fa", # key can also be "fasta_file" instead of "fasta"
#   "n": 10,
#   "delete_new_mod_bam_made_on_the_fly": false,
#   "genome_size_bp_optional": 0,
#   "no_jobs_for_head_on_and_co_directional": true,
#   "no_jobs_for_all": true,
#   "no_jobs_for_head_on": true,
#   "no_jobs_for_co_directional": true,
#   "pause_sensitivity_bedgraphs_prefix_optional": "/path/to/file",
#   "restrict_fork_direction_optional": "L",
#   "random_string_job_name": "job_name",
#   "regions_of_interest": [
#     {
#       "feature": "tRNA",
#       "bed_file": "/path/to/tRNA.bed",
#       "division": "NA"
#     },
#     {
#       "feature": "tRNA",
#       "bed_file": "/path/to/tRNA.1.bed",
#       "division": "01"
#     },
#     {
#       "feature": "tRNA",
#       "bed_file": "/path/to/tRNA.2.bed",
#       "division": "02"
#     },
#     {
#       "feature": "gene",
#       "bed_file": "/path/to/gene.bed",
#       "division": "NA"
#     },
#     {
#       "feature": "gene",
#       "bed_file": "/path/to/gene.1.bed",
#       "division": "01"
#     },
#     {
#       "feature": "gene",
#       "bed_file": "/path/to/gene.2.bed",
#       "division": "02"
#     }
#   ]
# }
# explanation of fields:
# - mod_bam: path to the modified bam file
# - forksense_dir: path to the forksense directory which has forks called by forksense.
# - pause_file: path to the pause file which contains all forks (whether they contain pauses or not) in a region
#               larger than the one being queried here (e.g.: the whole genome). this is a tab-separated file.
#               we are not describing the format of this file here, refer to scripts like convert_pause_file_to_bed.py.
# - op_dir: path to the output directory where the output files will be written. this is preferably new/empty as we are
#           going to write many files with some naming convention and it's better to not have any conflicts.
# - dataset: name of the dataset e.g.: 20180927_LAM_ONT_SC_wt60_a17aaa2
# - analysis_label: label for the analysis e.g.: 12jun23. the user can use any convenient label here; it's mostly
#                   for their own reference.
# - mod_bam_left: path to the modified bam file for left fork probabilities called by forkSense
# - mod_bam_right: path to the modified bam file for right fork probabilities called by forkSense
# - fasta: path to the fasta file containing the reference genome. For compatibility with other scripts, you can also
#          use the key "fasta_file" instead of "fasta".
# - n: number of typical reads to plot for each bed file. this is optional and defaults to 0 i.e. no reads are plotted.
#      NOTE: since we are going over many bed files and splitting by co-directional etc. and we want n reads for each
#            of these, the total number of reads plotted will be n * number of bed files * number of splits, which
#            could easily be a large number, so be careful and do not start a large number of plotting jobs.
# - delete_new_mod_bam_made_on_the_fly: (optional, default false) if true, delete the new mod_bam file made on the fly
#                                       This is the new mod bam file formed by intersecting the input mod bam with the
#                                       input bed file. If you are making many pause reports, you will make many new
#                                       mod bam files. If you don't delete them, they will take up a lot of space.
#                                       If you are only making one pause report, you can set this to false and keep
#                                       the new mod bam file.
# - genome_size_bp_optional: (optional, default 0) genome size in base pairs. Remember to include as many rDNAs as you
#                            think necessary and exclude mitochondrial genome if need be. Due to these assumptions, we
#                            do not calculate this automatically from the fasta file.
#                            By default or if set to 0, program will assume this information is not provided.
# - no_jobs_for_head_on_and_co_directional: (optional, default false) if true, do not get two additional pause reports
#                                           by splitting for head-on and co-directional conflicts.
# - no_jobs_for_all: (optional, default false) if true, do not get a pause report for all forks.
# - no_jobs_for_head_on: (optional, default false) if true, do not get a pause report for head-on forks.
# - no_jobs_for_co_directional: (optional, default false) if true, do not get a pause report for co-directional forks.
#                                The no_jobs* parameters all behave similarly.
# - pause_sensitivity_bedgraphs_prefix_optional: (default "" i.e. unused) prefix for pause sensitivity bedgraphs.
#                                               See bed_region_pause_report.sh for explanation.
#                                               See run_get_sgm_pause_sensitivities.sh for how to generate such files
#                                               or other details. By default, this field is not used and any associated
#                                               calculations are not performed.
# - restrict_fork_direction_optional: (default "" i.e. no restriction) restrict fork direction to "L"/"R"/"lead"/"lag"
#                                     i.e. only use left, right, leading or lagging forks respectively in any
#                                     calculation. We know that in reality forks perform both leading and lagging strand
#                                     synthesis. As our data is _single-molecule_, we _can_ classify forks as leading
#                                     or lagging. Default is no restriction.
# - random_string_job_name: (default "bedRegPauseReport") a random string to be used as the job name for the slurm jobs
#                           of each bed_region_pause_report.sh job.
# - regions_of_interest: see the explanation below.

# Region of interest specification
# --------------------------------
# In the regions_of_interest field, the user specifies the feature, the bed file(s) and the division in a list.
# Only the bed file is used by the program for calculations. The feature and the division are for labeling purposes
# only, so the user can specify any suitable string for these.
# An example entry is:
# {
#   "feature": "tRNA",
#   "bed_file": "/path/to/tRNA.1.bed",
#   "division": "01"
# }
# The interpretation here is that tRNAs have been divided up into many parts across several bed files and the user
# is specifying only the first division here and wants to get a pause report for it.
# Another example entry is:
# {
#   "feature": "tRNA",
#   "bed_file": "/path/to/tRNA.bed",
#   "division": "NA"
# }
# The interpretation here is that the bed file contains all the tRNAs without any division and we want a pause report
# for all the tRNAs.
# NOTE: do not use the same feature-division combination more than once in the list, as this will result in overwriting
#       of files.

# Overall logic
# -------------
# For each entry in the regions_of_interest list, we do three runs of bed_region_pause_report.sh, once for all
# forks, once for co-directional forks and once for head-on forks.
# So, if you have N entries in the regions_of_interest list, you will run 3*N jobs.

# outputs
# -------
# Multiple bed_region_pause_report.sh jobs are run.
# Refer to that script for details on the output files.
# Some messages are written to standard output/SBATCH output log file for logging purposes.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -jq

# load configuration
source config.sh

# print node information
bash print_node_information.sh

# assign arguments to variables
input_json=${1:-}

# check that the correct number of arguments were provided
if [ "$#" -lt 1 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash $0 input_file.json"
    >&2 echo "For more details on what these parameters mean, see the comments in the script."
    exit 1;
fi

# check that the input json file exists
if [ ! -f "$input_json" ]; then
    >&2 echo "ERROR: input json file $input_json does not exist"
    exit 1;
fi

# get the number of regions of interest
len_roi=$(< "$input_json" jq '.regions_of_interest | length')

# if there are no regions of interest, exit
if [ "$len_roi" -eq 0 ]; then
    >&2 echo "ERROR: no regions of interest specified"
    >&2 echo "Please specify at least one region of interest in the input json file"
    exit 1;
fi

# extract the op_dir from the input json file and make a temp directory in it
op_dir=$(< "$input_json" jq -r '.op_dir')
tmp_dir="$op_dir"/temp_"$(openssl rand -hex 6)"_made_on_"$(date +%Y%m%d)"_delete_later
mkdir -p "$tmp_dir"

# get the random string job name
random_string_job_name=$(< "$input_json" jq -r '.random_string_job_name // "bedRegPauseReport"')

# loop through the regions of interest and create a list of bed files, feature and division
bed_files=()
features=()
divisions=()
for counter in $(seq 1 "$len_roi")
do
  bed_files+=("$(< "$input_json" jq -r '.regions_of_interest['"$counter-1"'].bed_file')")
  features+=("$(< "$input_json" jq -r '.regions_of_interest['"$counter-1"'].feature')")
  divisions+=("$(< "$input_json" jq -r '.regions_of_interest['"$counter-1"'].division')")
done

# check that the bed files are unique
if [ "$(printf '%s\n' "${bed_files[@]}" | sort | uniq -d | wc -l)" -ne 0 ]; then
    >&2 echo "ERROR: duplicate bed files specified"
    >&2 echo "Please specify unique bed files in the input json file"
    exit 1;
fi

# check that the bed files exist
for bed_file in "${bed_files[@]}"
do
  if [ "$bed_file" == "null" ]; then
    >&2 echo "ERROR: bed file not specified in one of the regions of interest"
    exit 1;
  fi
  if [ ! -f "$bed_file" ]; then
    >&2 echo "ERROR: bed file $bed_file does not exist"
    exit 1;
  fi
done

# check that the features and divisions list do not contain null
for feature in "${features[@]}"
do
  if [ "$feature" == "null" ]; then
    >&2 echo "ERROR: feature not specified in one of the regions of interest"
    exit 1;
  fi
done

for division in "${divisions[@]}"
do
  if [ "$division" == "null" ]; then
    >&2 echo "ERROR: division not specified in one of the regions of interest"
    exit 1;
  fi
done

# get a list where each entry is the concatenation of each entry in the feature and division list
concat_feature_division=()
for counter in $(seq 1 "$len_roi")
do
  concat_feature_division+=("${features[$counter-1]}${divisions[$counter-1]}")
done

# check that the feature + division concatenations are unique
if [ "$(printf '%s\n' "${concat_feature_division[@]}" | sort | uniq -d | wc -l)" -ne 0 ]; then
    >&2 echo "ERROR: duplicate feature + division combinations specified"
    >&2 echo "Please specify unique feature + division combinations in the input json file"
    exit 1;
fi

# create a file that will contain all submitted jobs
submitted_jobs=$(mktemp --tmpdir="$tmp_dir" submitted_jobs.XXXXXX)

# Define a launch job record function that records the job launch in the submitted_jobs file and prints
# some information to the standard output
launch_job_record() {
    modified_json=${1:-}
    index=${2:-NA}

    if [ ! -f "$modified_json" ]; then
        >&2 echo "ERROR: modified json file $modified_json does not exist"
        exit 1;
    fi

    # Append the command to submitted_jobs file
    echo "bash bed_region_pause_report.sh $modified_json" >> "$submitted_jobs"

    # Print some information to the standard output
    echo "[LOG] Launched job for the following json object $modified_json in position $index in job array:"
    cat "$modified_json"
}

# loop through the regions of interest
job_counter=1
for counter in $(seq 1 "$len_roi")
do

  # separate out each bed file into its own json file
  < "$input_json" jq 'del(.regions_of_interest) + .regions_of_interest['"$counter-1"']' |\
    jq '. + {"prefix_plot_option": true}' > "$tmp_dir"/temp"$counter".json

  # retrieve all the no jobs flags
  no_jobs_for_head_on_and_co_directional=$(< "$input_json" jq -r '.no_jobs_for_head_on_and_co_directional // false')
  no_jobs_for_all=$(< "$input_json" jq -r '.no_jobs_for_all // false')
  no_jobs_for_head_on=$(< "$input_json" jq -r '.no_jobs_for_head_on // false')
  no_jobs_for_co_directional=$(< "$input_json" jq -r '.no_jobs_for_co_directional // false')

  # first, do a job for all forks
  if [ ! "$no_jobs_for_all" == "true" ]; then
    prefix_str="${features[$counter-1]}_${divisions[$counter-1]}_all"
    modified_json="$tmp_dir"/temp"$counter"_all.json
    jq '. + {"relative_direction_option": "all", "prefix": "'"$prefix_str"'"}' "$tmp_dir"/temp"$counter".json \
      > "$modified_json"
    launch_job_record "$modified_json" "$job_counter"
    job_counter=$((job_counter+1))
  else
    echo "[LOG] Skipping job for all forks i.e. indiscriminate of direction"
  fi

  # if no jobs for head-on and co-directional, skip the next steps
  if [ "$no_jobs_for_head_on_and_co_directional" == "true" ]; then
    continue
  else
    echo "[LOG] Skipping job for head-on and co-directional forks"
  fi

  # run jobs to get separate pause reports
  if [ ! "$no_jobs_for_head_on" == "true" ]; then
    prefix_str="${features[$counter-1]}_${divisions[$counter-1]}_ho"
    modified_json="$tmp_dir"/temp"$counter"_ho.json
    jq '. + {"relative_direction_option": "head-on", "prefix": "'"$prefix_str"'"}' "$tmp_dir"/temp"$counter".json \
      > "$modified_json"
    launch_job_record "$modified_json" "$job_counter"
    job_counter=$((job_counter+1))
  else
    echo "[LOG] Skipping job for head-on forks"
  fi

  if [ ! "$no_jobs_for_co_directional" == "true" ]; then
    prefix_str="${features[$counter-1]}_${divisions[$counter-1]}_cd"
    modified_json="$tmp_dir"/temp"$counter"_cd.json
    jq '. + {"relative_direction_option": "co-directional", "prefix": "'"$prefix_str"'"}' "$tmp_dir"/temp"$counter".json \
      > "$modified_json"
    launch_job_record "$modified_json" "$job_counter"
    job_counter=$((job_counter+1))
  else
    echo "[LOG] Skipping job for co-directional forks"
  fi

done

# submit all jobs
n_jobs=$((job_counter-1))
if [ "$n_jobs" -eq 0 ]; then
  >&2 echo "ERROR: no jobs to submit"
  exit 1;
else
  echo "[LOG] Submitting $n_jobs jobs with the array name $random_string_job_name"
  sbatch --mem=10G -c 1 --time=3:59:59 -J "$random_string_job_name" --mail-user= -p ei-short  \
  --array=1-"$n_jobs" job_array_converter.sh "$submitted_jobs"
fi