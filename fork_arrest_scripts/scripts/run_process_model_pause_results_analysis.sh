#!/bin/bash

#SBATCH --mem=100G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J processPauseAnalysis
#SBATCH --mail-type=END,FAIL
#SBATCH --time 167:59:59
#SBATCH --constraint=""

# goal
# -----
# run the process model pause results analysis script.

# usage
# -----
# sbatch run_process_model_pause_results_analysis.sh input_file output_file modbam_file modbam_left modbam_right
# input_file: input file
# output_file: output file
# modbam_file: modbam file containing analogue modification probabilities per thymidine per read
# modbam_left: modbam file containing fork probabilities per thymidine per read for left forks
# modbam_right: modbam file containing fork probabilities per thymidine per read for right forks

# preamble
# --------

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load python
source load_package.sh -python

# get configuration
source config.sh

# check that only five inputs are present on the command line
if [ $# -ne 5 ]; then
  echo "usage: sbatch run_process_model_pause_results_analysis.sh input_file output_file modbam modbam_left modbam_right"
  exit 1
fi

# set input and output files
ip_file=$1
op_file=$2
modbam=$3
modbam_left=$4
modbam_right=$5

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

# set the script information
infoStr="# from commit ${COMMITSTR} generated at ${TIMENOW}"

{
  # output script information
  echo "$infoStr"

  # perform analysis of pause results
  # shellcheck disable=SC2016
  < "$ip_file" python process_model_pause_results_analysis.py\
    --modbam "$modbam"\
    --modbam_left "$modbam_left" --modbam_right "$modbam_right"\
    --keep_cols_requested keep_exp_gof_v2,keep_fl,keep_fp,keep_ep,keep_lp,keep_fs
} > "$op_file"
