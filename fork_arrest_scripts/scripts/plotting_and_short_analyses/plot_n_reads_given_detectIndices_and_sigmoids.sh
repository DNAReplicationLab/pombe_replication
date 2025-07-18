#!/bin/bash

#SBATCH --mem=10G
#SBATCH -c 3
#SBATCH -p ei-short
#SBATCH -J plot_n_reads_given_detectIndices_and_sigmoids
#SBATCH --mail-type=NONE
#SBATCH --mail-user=
#SBATCH --time 0:39:59

# introduction
# --------------
# We want to plot the raw and windowed data of analogue probabilities, along with model fit information
# of a fork. The model we are using here is a sigmoid fit to the raw, non-windowed data where the high and
# low levels are fixed, but the width and the location of the center point of the sigmoid are fitting parameters.
# Additionally, we want to show the rest of the read and the features on it, such as forks, pauses, origins,
# terminations.

# usage
# --------
# bash <program_name.sh> $detect_index_sigmoid_file $mod_bam $forksense_dir $fasta_file $mod_bam_left $mod_bam_right \
#   $op_dir [$window_size] [$threshold]
# NOTE: The "\" symbol means the command is continued on the next line.
# NOTE: <> means substitute whatever is within the angle brackets suitably.
# NOTE: [] means the argument is optional. If you want to specify the argument, remove the brackets and
#       write the argument. If you don't want to specify the argument, just don't write anything.
#       To specify an optional argument, you must specify all arguments to the left of it as well.
# detect_index_sigmoid_file:
#   space-separated file with no headers. The first column is the detect index,
#   the second (and possibly subsequent) column is a string of a particular format containing
#   the best-fit sigmoid's parameters.
#   * detect index means readID_contig_start_end_orientation_direction_startFork_endFork.
#     readID is the read id, contig is the contig name, start and end are the start and end of the alignment to the
#     reference genome, orientation is the orientation of the read (fwd/rev), direction is the direction of the fork
#     (L/R), startFork and endFork are the start and end of the fork along reference coordinates,
#     with startFork < endFork irrespective of the fork direction. All units are bp.
#   * The string of parameters is of the format: high_low_offset_width_start_end.
#     high and low are the high and low levels of the sigmoid,
#     offset is the center point of the sigmoid given by coordinates on the reference in bp,
#     width is the width of the sigmoid (positive if fork direction is R, negative if fork direction is L),
#     start and end are where to plot the sigmoid from and to on the reference.
#     All units are bp.
#    * There can be many columns of parameter strings, and the number of parameter columns per row can be variable.
# mod_bam: the .mod.bam file that contains analogue modification probabilities of the reads of interest
# forksense_dir: directory with forksense files that contain origin, left fork, right fork, termination information
#                in the standard forksense format and file-naming convention.
# fasta_file: fasta file of the reference genome
# mod_bam_left: the .mod.bam file that contains left-fork probabilities of the reads of interest
# mod_bam_right: the .mod.bam file that contains right-fork probabilities of the reads of interest
# op_dir: output directory where all plots/data are sent.
# window_size: (optional, default 300) window size in thymidines for windowing analogue probabilities.
# threshold: (optional, default 0.5) threshold above which thymidines are called as BrdU.

# stop execution if any command fails
set -e

if [ "$#" -lt 7 ]
then
  echo "Usage: bash <program_name.sh> detect_index_sigmoid_file mod_bam forksense_dir fasta_file mod_bam_left mod_bam_right op_dir [window_size] [threshold]"
  exit
fi

# set input directories, files
input_file=$1
mod_bam=$2
forksense_dir=$3
fasta_file=$4
mod_bam_left=$5
mod_bam_right=$6

# set output directory
op_dir=$7

# set window size
window_size=${8:-300}

# set threshold
threshold=${9:-0.5}

# load samtools
pwd=$(pwd)
cd ..
source load_package.sh -samtools
cd "$pwd"

# check if input file exists
if [ ! -f "$input_file" ]
then
  echo "Input file $input_file does not exist. Exiting."
  exit
fi

# check if mod_bam file exists
if [ ! -f "$mod_bam" ]
then
  echo "Input file $mod_bam does not exist. Exiting."
  exit
fi

# check if forksense_dir exists
if [ ! -d "$forksense_dir" ]
then
  echo "Input directory $forksense_dir does not exist. Exiting."
  exit
fi

# check if fasta_file exists
if [ ! -f "$fasta_file" ]
then
  echo "Input file $fasta_file does not exist. Exiting."
  exit
fi

# mod_bam_left and mod_bam_right are optional. So we are not checking if they exist here.

# make output directory if it does not exist
mkdir -p "$op_dir"
op_dir=$(cd "$op_dir" && pwd)

# extract read ids of interest
awk -v IFS="_" '{print $1}' "$input_file" > "${op_dir}/reads.txt"

# check that all reads are unique; it is possible for read ids to repeat in the input file.
# We cannot deal with this situation right now, so we print an error and exit.
n_repeated_reads=$(sort "${op_dir}/reads.txt" | uniq -d | wc -l)
if [ "$n_repeated_reads" -ne 0 ]
then
  echo "Read ids in the input file are not unique. Exiting."
  exit 1;
fi

# Function to create subsets of bam files
create_bam_subset() {
  local bam_file="$1"
  local subset_prefix="$2"
  local bam_subset="/dev/null"

  if [[ -f "$bam_file" ]]; then
    bam_subset="${op_dir}/${subset_prefix}.bam"
    samtools view -b -N "${op_dir}/reads.txt" "$bam_file" > "${bam_subset}.tmp"
    samtools sort -o "$bam_subset" "${bam_subset}.tmp"
    samtools index "$bam_subset"
    rm "${bam_subset}.tmp"
  fi

  echo "$bam_subset"
}

# * Make subsets of bam files wherever possible
# ** Create temporary files to capture the output of background tasks
tmp_file_1=$(mktemp -p "$op_dir")
tmp_file_2=$(mktemp -p "$op_dir")
tmp_file_3=$(mktemp -p "$op_dir")

# ** Make subsets of bam files wherever possible and capture their output in temporary files
(create_bam_subset "$mod_bam" "mod_bam_subset" > "$tmp_file_1") &
(create_bam_subset "$mod_bam_left" "mod_bam_left_subset" > "$tmp_file_2") &
(create_bam_subset "$mod_bam_right" "mod_bam_right_subset" > "$tmp_file_3") &

# ** Wait for all background tasks to finish
wait

# ** Read the outputs from the temporary files back into variables
mod_bam_subset=$(<"$tmp_file_1")
mod_bam_left_subset=$(<"$tmp_file_2")
mod_bam_right_subset=$(<"$tmp_file_3")

# ** Remove the temporary files
rm "$tmp_file_1" "$tmp_file_2" "$tmp_file_3"

# get a randomly generated string, will be used for job names
job_name=plot"$(openssl rand -hex 6)"

# read data line by line from file and start plotting jobs
while IFS=' ' read -r -a line_data; do

  IFS='_' read -r -a detect_index_data <<< "${line_data[0]}"

  # extract data from the detect index
  read_id="${detect_index_data[0]}"
  contig="${detect_index_data[1]}"
  start="${detect_index_data[2]}"
  end="${detect_index_data[3]}"
  start_fork="${detect_index_data[6]}"
  end_fork="${detect_index_data[7]}"
  center_point="$(echo "scale=0; (${start_fork} + ${end_fork}) / 2" | bc)"

  # get feature information to temporary file
  {
    echo "# coordinates on the reference genome: $contig $start $end"
    echo "# detectIndex: ${line_data[0]}"
    echo "# sigmoid parameter(s): ${line_data[*]:1}"
    bash get_feature_information.sh "$read_id" /dev/null "$forksense_dir"
  } > "$op_dir"/features_"$read_id"

  # prepare parameter string
  parameter_string=""
  loop_counter=0

  # loop through the line_data variable after the first entry
  for sig_param in "${line_data[@]:1}"; do
    # if loop counter is zero, then just assign the sigmoid parameters to the parameter string
    if [ "$loop_counter" -eq 0 ]
    then
      parameter_string="$sig_param "
    else
      # else add the sigmoid parameters with a break in between to show these are separate curves
      parameter_string+="${center_point}_NA $sig_param "
    fi
    ((++loop_counter))
  done

  # remove the last character (the extra space) from the parameter string
  parameter_string=${parameter_string%?}

  # plot analogue data along the fork with sigmoid fit
  # windowing is done such that the center of the fork falls on a window boundary
  # shellcheck disable=SC2086
  sbatch --mail-user= -J "$job_name" plot_read.sh "$mod_bam_subset" "$read_id" "$contig" "$start" "$end" "$op_dir"\
            "$window_size" "$threshold" "$center_point" "$fasta_file" $parameter_string\
            "$op_dir"/features_"$read_id"

  # plot forksense left
  if [[ -f "$mod_bam_left_subset"  ]]; then
    sbatch --mail-user= -J "$job_name" limited_plot_read.sh "$mod_bam_left_subset" "$read_id" 0 "$op_dir" \
      forkSense_left;
  fi

  # plot forksense right
  if [[ -f "$mod_bam_right_subset"  ]]; then
    sbatch --mail-user= -J "$job_name" limited_plot_read.sh "$mod_bam_right_subset" "$read_id" 0 "$op_dir" \
      forkSense_right;
  fi

done < "$input_file"

# collate plots into a pdf
sbatch -p ei-medium --time=1:30:00 --mail-user= --dependency=singleton --job-name="$job_name"\
    convert_plots_to_latex.sh /dev/null "$op_dir" "$forksense_dir" "$op_dir" forkSense features