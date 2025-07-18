#!/bin/bash

#SBATCH --mem=100G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J parallelProcessPauseAnalysis
#SBATCH --mail-type=END,FAIL
#SBATCH --time 71:59:59
#SBATCH --constraint=""

# goal
# -----
# break input pause file into pieces and run the process model pause results analysis script on each piece,
# and stitch them back together.

# usage
# -----
# sbatch run_parallel_run_process_model_pause_results_analysis.sh input_file output_file modbam_file modbam_left\
#   modbam_right
# NOTE: \ at end of line is for line continuation and need not be part of the command if it is all on one line.
# input_file: input file
# output_file: output file
# modbam_file: modbam file containing analogue modification probabilities per thymidine per read
# modbam_left: modbam file containing fork probabilities per thymidine per read for left forks
# modbam_right: modbam file containing fork probabilities per thymidine per read for right forks

# stop execution if any command fails
set -e

# load packages
source load_package.sh -miller

# check that only five inputs are present on the command line
if [ $# -ne 5 ]; then
  echo "usage: sbatch run_parallel_run_process_model_pause_results_analysis.sh input_file output_file modbam modbam_left modbam_right"
  exit 1
fi

# set input and output files
ip_file=$1
op_file=$2
modbam=$3
modbam_left=$4
modbam_right=$5
n_proc=50 # number of independent slurm processes to run

# check that input files exist
if [ ! -f "$ip_file" ]; then
  echo "input file $ip_file does not exist"
  exit 1
fi

# check that bam files exist
if [ ! -f "$modbam" ]; then
  echo "modbam file $modbam does not exist"
  exit 1
fi

if [ ! -f "$modbam_left" ]; then
  echo "modbam file $modbam_left does not exist"
  exit 1
fi

if [ ! -f "$modbam_right" ]; then
  echo "modbam file $modbam_right does not exist"
  exit 1
fi

# check that bam files are indexed
if [ ! -f "${modbam}.bai" ]; then
  echo "modbam file $modbam is not indexed"
  exit 1
fi

if [ ! -f "${modbam_left}.bai" ]; then
  echo "modbam file $modbam_left is not indexed"
  exit 1
fi

if [ ! -f "${modbam_right}.bai" ]; then
  echo "modbam file $modbam_right is not indexed"
  exit 1
fi

# get a temporary directory
tmp_dir="$op_file".temp_"$(openssl rand -hex 6)"_made_on_"$(date +%Y%m%d)"_delete_later
mkdir -p "$tmp_dir"

# split the input file into many pieces
mlr --tsv --skip-comments split -m "$n_proc" --prefix "$tmp_dir"/input_split "$ip_file"

# create a file that will contain all submitted jobs
submitted_jobs=$(mktemp --tmpdir="$tmp_dir" submitted_jobs.XXXXXX)

# Define the sbatch_record function that stores the command to be submitted in a file
sbatch_record() {
    # Append the command to submitted_jobs file
    echo "bash $*" >> "$submitted_jobs"
}

# perform analysis of pause results
# each job is a slurm process

for count in $(seq 1 1 $n_proc); do

  ip_file_piece="$tmp_dir"/input_split_"$count".tsv

  if [ ! -f "$ip_file_piece" ]; then
    echo "input file piece $ip_file_piece does not exist"
    exit 1
  fi

  sbatch_record run_process_model_pause_results_analysis.sh\
    "$ip_file_piece" "${ip_file_piece}.processed" "$modbam" "$modbam_left" "$modbam_right"

done

# submit all jobs
n_jobs=$(wc -l < "$submitted_jobs")
launchJob=$(sbatch --mem=10G -c 1 --time=11:59:59 -J processPauseAnalysis --mail-user= -p ei-short  \
  --array=1-"$n_jobs" job_array_converter.sh "$submitted_jobs")
jid=${launchJob##* }

# wait for all slurm jobs to finish
# we set time to 30 minutes although the job is sleep 10 (seconds) because
# some times even slurm jobs with simple commands like sleep 10 take unusually long to run.
sbatch -p ei-short --wait --time=30:00 --mail-user= --dependency=afterok:"$jid" --job-name="sleep_10" --wrap="sleep 10"

{

  # output any comments in the input file
  head -n 100 "$ip_file" | grep '^#'

  # output any comments in each file.
  for file in "$tmp_dir"/input_split_[0-9]*.tsv.processed; do
    head -n 100 "$file" | grep '^#'
  done

  # join the pieces
  mlr --tsv --skip-comments cat "$tmp_dir"/input_split_[0-9]*.tsv.processed

} > "$op_file"

# remove the temporary directory
rm -rf "$tmp_dir"
