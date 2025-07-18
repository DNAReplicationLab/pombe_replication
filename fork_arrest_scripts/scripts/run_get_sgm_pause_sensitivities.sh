#!/bin/bash

#SBATCH --mem-per-cpu=100G
#SBATCH -c 3
#SBATCH -p ei-medium
#SBATCH -J runGetSgmPauseSensitivities
#SBATCH --mail-type=END,FAIL
#SBATCH --time=47:59:59
#SBATCH --constraint=""

# goal
# -----
# Use detected pauses to obtain likelihood of detecting pauses across the genome.
# For logic used, refer to the individual scripts used below.

# usage and inputs
#-----------------
# sbatch run_get_sgm_pause_sensitivities.sh [-w] pause_file keep_columns fasta_file output_prefix window_size_optional\
#   decimate_T_optional
# use sbatch. do not use bash unless you are on a node with enough memory, time etc.
# NOTE: \ means that the command continues on the next line.
# NOTE: [] means that the argument is optional. If you don't want to set it, just remove all text between the brackets
#       including the brackets. If you want to set it, remove the brackets and write the text within.
# -w: optional (default: overwrite files). If set, existing output files will not be overwritten.
#     This is useful if you have already obtained the single base bed file for example, and want to regenerate
#     the windowed sensitivity files with a different window size (w stands for the write in overwrite, not window
#     size).
# pause_file: pause file in our standard tab-separated pause format.
# keep_columns: comma-separated list of columns starting with "keep" to be used while assigning pause probabilities
# fasta_file: reference genome in fasta format. Must have a corresponding index file called "$fasta_file".fai
# output_prefix: prefix for output files. so if you want to name your output file "my_output_file.bedgraph",
#                and put it in the folder "/my/output/folder", then you should set output_prefix to
#                "/my/output/folder/my_output_file".
# window_size_optional: window size in bp to use for calculating pause probabilities. Default is 1000.
# decimate_T_optional: only output probability vs fork coordinate every decimate_T_optional thymidines.
#                      default is 20. saves compute time.
# NOTE: if you want further explanation of the parameters, see the comments in the scripts called below.
#       We've intentionally kept the description here brief as it is a wrapper script.

# outputs
# -------
# - A bed file and seven bedgraph files beginning with the output_prefix are generated.
# - Bed file contains per-read per-base pause sensitivity, and the seven bedgraphs are summations
#   of the sensitivity along the genome across all, left (L), right (R), L+, L-, R+, R- forks.
# - Some log messages are printed to stdout.
# - As this is a wrapper script that exists just to check inputs and calls other scripts,
#   we have not included a detailed description of the output files, which are in the respective scripts.
# - IMPORTANT: In the all fork file, the sensitivity column is normalized to add to the total number of pauses in the
#              pause file. In the left and right fork file, the same normalization factor as the all fork file is used.
#              This means the sum of the sensitivity column in the left file is NOT the total number of pauses detected
#              on left forks and the same is true for the right fork file. The same normalization value ensures
#              we can add intervals easily across different files.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load configuration
source config.sh

# print node information
bash print_node_information.sh

# function to set calling script information
insert_calling_script_header() {
  sed '1i'\
'# from commit '"${COMMITSTR:-NA}"' generated at '"${TIMENOW:-NA}"' by '"${config[name]:-NA}"' <'"${config[email]:-NA}"'>\n'\
'# script: '"$0"'\n'\
'# arguments: '"$*"'\n'\
"# slurm job name: ${SLURM_JOB_NAME:-NA}"
}

# process flags
is_do_not_overwrite=0

while getopts ":w" opt; do
  case $opt in
    w)
      is_do_not_overwrite=1
      echo "INFO: existing output files will not be overwritten"
      ;;
    \?)
      >&2 echo "ERROR: invalid option: $OPTARG"
      exit 1
      ;;
  esac
done

# shift arguments to exclude parsed options
shift $((OPTIND-1))

# assign arguments to variables
pause_file=${1:-}
keep_columns=${2:-}
fasta_file=${3:-}
output_prefix=${4:-}
window_size=${5:-1000}
decimate_T=${6:-20}

# check that the correct number of arguments were provided
if [ "$#" -lt 4 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: sbatch run_get_sgm_pause_sensitivities.sh pause_file keep_columns fasta_file output_prefix \ "
    >&2 echo "              window_size_optional decimate_T_optional"
    >&2 echo "use sbatch. do not use bash unless you are on a node with enough memory, time etc."
    >&2 echo "NOTE: \ means that the command continues on the next line."
    >&2 echo "pause_file: pause file in our standard tab-separated pause format."
    >&2 echo "keep_columns: comma-separated list of columns starting with \"keep\" to be used while assigning pause probabilities"
    >&2 echo "fasta_file: reference genome in fasta format. Must have a corresponding index file called \"\$fasta_file\".fai"
    >&2 echo "output_prefix: prefix for output files. see comments in script for more details on this parameter."
    >&2 echo "window_size_optional: window size in bp to use for calculating pause probabilities. Default is 1000."
    >&2 echo "decimate_T_optional: only output probability vs fork coordinate every decimate_T_optional thymidines."
    >&2 echo "                     default is 20. saves compute time."
    >&2 echo "The suffix _optional means that the parameter is optional."
    exit 1;
fi

# check that window size is an integer
if ! [[ "$window_size" =~ ^[0-9]+$ ]]; then
    >&2 echo "ERROR: window size must be an integer"
    exit 1;
fi

# check that window size is above 0
if [ "$window_size" -le 0 ]; then
    >&2 echo "ERROR: window size must be above 0"
    exit 1;
fi

# check if fasta index file exists
if [ ! -f "$fasta_file".fai ]; then
    >&2 echo "ERROR: fasta index file $fasta_file.fai does not exist"
    exit 1;
fi

# check that pause file exists
if [ ! -f "$pause_file" ]; then
    >&2 echo "ERROR: pause file $pause_file does not exist"
    exit 1;
fi

# check that decimate_T is an integer
if ! [[ "$decimate_T" =~ ^[0-9]+$ ]]; then
    >&2 echo "ERROR: decimate_T must be an integer"
    exit 1;
fi

# check that decimate_T is above 0
if [ "$decimate_T" -le 0 ]; then
    >&2 echo "ERROR: decimate_T must be above 0"
    exit 1;
fi

# set filename for intermediate steps
probabilistic_single_base_file="$output_prefix".pause_prob_per_read_per_single_base.bed

# perform steps to get pause sensitivity
if [ ! -f "$probabilistic_single_base_file" ] || [ "$is_do_not_overwrite" -eq 0 ]; then
  echo "INFO: converting pause file to probabilistic single base bed file"
  bash convert_pause_file_to_probabilistic_single_base_bed.sh "$pause_file"\
    "$keep_columns" "$fasta_file" "$decimate_T" > "$probabilistic_single_base_file"
else
  echo "INFO: probabilistic single base bed file already exists. Skipping conversion"
fi

fork_direction_regex=("^L$" "^R$" "^L|R$");
fork_direction_suffix=("left" "right" "all");

for count in {0..2}; do
  current_output_file="$output_prefix".pause_sensitivity."${fork_direction_suffix[$count]}".bedgraph
  if [ -f "$current_output_file" ] && [ "$is_do_not_overwrite" -eq 1 ]; then
    echo "INFO: Will not overwrite the pre-existing $current_output_file as -w flag was set."
    continue
  else
    echo "INFO: Ongoing creation of $current_output_file ..."
  fi
  {

    echo "# filter to retain ${fork_direction_suffix[$count]} forks"

    # print comments
    head -n 100 "$probabilistic_single_base_file" | grep '^#'

    # perform calculations
    grep -v '^#' "$probabilistic_single_base_file" |\
      awk -v r="${fork_direction_regex[$count]}" 'BEGIN{OFS="\t"}{if($8 ~ r){print $1, $2, $3, $4, $7, $6}}' |\
        bash sum_per_read_per_base_bed_to_genome_window_bed.sh - "$fasta_file".fai "$window_size"

  } | insert_calling_script_header "$@" > "$current_output_file" &
done

# get L+, L-, R+, R- pause sensitivities - these will be useful for analysis
# separating by leading and lagging strand synthesis

fork_direction_regex=("^L$" "^R$" "^L$" "^R$");
strand=("+" "+" "-" "-");
file_suffix=("left.plus" "right.plus" "left.minus" "right.minus");

for count in {0..3}; do
  current_output_file="$output_prefix".pause_sensitivity."${file_suffix[$count]}".bedgraph
  if [ -f "$current_output_file" ] && [ "$is_do_not_overwrite" -eq 1 ]; then
    echo "INFO: Will not overwrite the pre-existing $current_output_file as -w flag was set."
    continue
  else
    echo "INFO: Ongoing creation of $current_output_file ..."
  fi
  {

    echo "# filter to retain ${file_suffix[$count]} forks"

    # print comments
    head -n 100 "$probabilistic_single_base_file" | grep '^#'

    # perform calculations
    grep -v '^#' "$probabilistic_single_base_file" |\
      awk -v r="${fork_direction_regex[$count]}" -v q="${strand[$count]}" 'BEGIN{OFS="\t"}{if($8 ~ r && $6 == q){print $1, $2, $3, $4, $7, $6}}' |\
        bash sum_per_read_per_base_bed_to_genome_window_bed.sh - "$fasta_file".fai "$window_size"

  } | insert_calling_script_header "$@" > "$current_output_file" &
done

wait;