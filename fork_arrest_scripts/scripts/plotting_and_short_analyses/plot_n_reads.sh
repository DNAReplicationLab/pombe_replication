#!/bin/bash

#SBATCH --mem=10G
#SBATCH -c 4
#SBATCH -p ei-medium
#SBATCH -J lPlotNReads
#SBATCH --mail-type=END,FAIL
#SBATCH --time 0:59:59

# introduction
# --------------

# Given a list of read ids, make three plots per read id.
# Plot 1: plot raw and windowed data of analogue probabilities, along with model fit information
#         and fork feature annotations
# Plots 2, 3: plot raw left and right fork probabilities as called by fork sense

# usage
# --------

# the script takes many inputs, so please read through the comments.

# bash <program_name.sh> [-c] $read_id_file $mod_bam $mod_bam_left $mod_bam_right $forksense_dir $pause_file \
#     $fasta_file $op_dir [$bam_file] [$n_reads] [$window_size] [$threshold]
# NOTE: The "\" symbol means the command is continued on the next line.
# NOTE: <> means substitute whatever is within the angle brackets suitably.
# NOTE: [] means the argument is optional. If you want to specify the argument, remove the brackets and
#       write the argument. If you don't want to specify the argument, just don't write anything.
#       To specify the argument, you must specify all arguments to the left of it as well.
# NOTE: even if some file options are not marked as optional in the documentation using [], they can be set to /dev/null
#       if you don't want to use them. /dev/null is a special file in linux systems that means 'empty file'.
#       Please see the section "variations" to see how to do this. Do not use /dev/null for numeric arguments.
# NOTE: can use sbatch instead of bash to run the script if you are on a cluster.

# Explanation of inputs
# option -c: If included, plot the reference sigmoid for all forks where pauses were not found.
#            This is useful to see if the sigmoid fits the data well.
#            Not done by default.
#            If you want to specify this option, remove the square brackets and write -c.
# read_id_file: This can be a one-column text file, a three-column text file or a bed file. See below for details.
#               The program automatically detects which format the input file is in and processes it accordingly.
#               -------------------------
#               Option 1: one column file
#               -------------------------
#               a one column text file with no header that contains one read id per line.
#               It is advisable to use a small number of read ids, say 10-50.
#               Use 100s or 1000s of ids only if you know what you are doing.
#               ---------------------------
#               Option 2: Three column file
#               ---------------------------
#               optionally, can use two more columns that show a start and an end, using a space as a column separator.
#               if the start and end of a fork feature of that read id coincide with these, then that feature
#               will be highlighted.
#               the goal here is to highlight a fork or origin of interest when the read id has multiple features.
#               -----------------
#               Option 3: bed file
#               -----------------
#               a bed file with at least 4 columns: contig, start, end, read_id.
#               BED format means tab-separated with no column names.
#               If a feature has a start and an end coinciding with the start and end of the read id,
#               then that feature is highlighted, like in option 2.
#               -----------------
#               Option 4: pause file
#               -----------------
#               a pause file. This is the format in which pause data is stored in our workflows.
#               This is a tab-separated file with headers and must contain at least two columns:
#               detectIndex and pauseSite.
#               Columns starting with "keep*" are used to filter out pauses.
#               For more details on the format, see convert_pause_file_to_bed.py
# mod_bam: the .mod.bam file that contains analogue modification probabilities of the reads of interest
# mod_bam_left: mod bam file with left fork probabilities per thymidine of the reads of interest
# mod_bam_right: same as above but with right fork probabilities
# forksense_dir: directory with forksense files that contain origin, left fork, right fork, termination information.
# pause_file: file with pause information. This is a tab separated file with headers.
#             the script uses the columns detectIndex, keep*, pauseSite, paramStrLeft, paramStrRight.
#             * means a wildcard i.e. any column that starts with 'keep'.
#             only rows with all keep* values = True are considered as pauses.
# fasta_file: reference genome
# op_dir: output directory where all plots are sent.
#         it is advisable to use directories without any plots in them.
#         ignore this rule only if you've read through the script and you know what you are doing.
# bam_file: optional. If bam file is given, script plots query coordinate vs reference coordinate for each read.
#          This is useful to see if the read is aligned properly. Default is to not plot this.
# n_reads: optional. If specified, reads are shuffled and so many reads are selected.
#          If n_reads is greater than the number of reads in the input file, all reads are plotted.
#          Default is 400. If you want to plot more than 400 reads, you must specify this argument.
#          Also, if you want to plot more than 400, ask yourself why do you want to plot so many reads?
#          It might be better to plot a subset of reads.
#          BTW, we have set the default to 400 to cover a range of possibilities. It does not mean that
#          the script expects to see 100s of reads in the input file. It is just a default value.
#          To not flood the cluster, plot a small number of reads, say 10-50.
# window_size: (optional) window size in thymidines to perform averaging for the plot. Default is 300.
# threshold: (optional) threshold for calling a thymidine as BrdU. Default is 0.5.

# variations
# ----------

## you are only interested in analogue modifications vs coordinate along reads.
## you are not interested in anything else like forks, models etc.
## use the following command

# bash <program_name.sh> $read_id_file $mod_bam /dev/null /dev/null /dev/null /dev/null /dev/null $op_dir
## as you can see, all unused inputs have been replaced with /dev/null, a special file in linux systems that
## means 'empty file'.


# stop execution if any command fails
set -e

# parse options
while getopts "c" opt; do
  case $opt in
    c)
      plot_option="-c"
      ;;
    \?)
      >&2 echo "Invalid option: -$OPTARG"
      exit 1
      ;;
  esac
done

# shift the options off the arguments
shift $((OPTIND-1))

if [ "$#" -lt 8 ]
then
  echo "Incorrect number of arguments! "
  exit
fi

# set some directories
current_dir=$(pwd)
main_dir=$(cd .. && pwd)

# load packages
cd "$main_dir"
source load_package.sh -python -samtools -jq
cd "$current_dir"

# set input directories, files
read_id_file=$1
mod_bam=$2
mod_bam_left=$3
mod_bam_right=$4
forksense_dir=/dev/null
if [[ -d "$5" ]]; then
  forksense_dir=$(cd "$5"; pwd)
fi
pause_file=$6
fa=$7

# check that the pause file is valid if it's specified
if [ -f "$pause_file" ]; then
  cd "$main_dir";
  if [ ! "$(< "$pause_file" python validate_pause_format.py)" == "valid" ]; then
    >&2 echo "Error: pause file is not valid."
    exit 1;
  fi
  cd "$current_dir";
fi

# set output directory, making it if it doesn't exist
mkdir -p "$8"
op_dir=$(cd "$8"; pwd)

# set bam file if it's given
bam_file=${9:-/dev/null}

# set the number of reads to plot if it's given
n_reads=${10:-400}

# set the window size if it's given
window_size=${11:-300}

# set the threshold if it's given
threshold=${12:-0.5}

# check if the read id file is a bed file
read_id_processed=$(mktemp -p "$op_dir")
cd "$main_dir"
{
  if [ "$(< "$read_id_file" python validate_bed_format.py --require-uuid --allow-float-score)" == "valid" ]; then
    < "$read_id_file" python print_valid_data_bed_lines.py | awk 'BEGIN{OFS=" "}{print $4, $2, $3}'
  elif [ "$(< "$read_id_file" python validate_pause_format.py)" == "valid" ]; then
    < "$read_id_file" python convert_pause_file_to_bed.py |\
       python print_valid_data_bed_lines.py | awk 'BEGIN{OFS=" "}{print $4, $2, $3}'
  else
    grep -v "^#" "$read_id_file"
  fi
} | shuf | head -n "$n_reads" > "$read_id_processed"
cd "$current_dir"

# check that all reads are unique; it is possible for read ids to repeat in the input file.
# We cannot deal with this situation right now, so we print an error and exit.
n_repeated_reads=$(sort "$read_id_processed" | uniq -d | wc -l)
if [ "$n_repeated_reads" -ne 0 ]
then
  echo "Read ids in the input file are not unique. Exiting."
  exit 1;
fi

echo "Plotting the following reads:"
cat "$read_id_processed"

# make subset of bam files
cd "$main_dir"
bash subset_four_bam_files.sh "$bam_file" "$mod_bam" "$mod_bam_left" "$mod_bam_right" "$read_id_processed" "$op_dir" \
  >  "$op_dir"/bam_subset_paths.tmp.json
cd "$current_dir"

# ** Read the outputs from the temporary files back into variables
bam_file_subset=$(jq -r '.bam_subset' "$op_dir"/bam_subset_paths.tmp.json)
mod_bam_subset=$(jq -r '.mod_bam_subset' "$op_dir"/bam_subset_paths.tmp.json)
mod_bam_left_subset=$(jq -r '.mod_bam_left_subset' "$op_dir"/bam_subset_paths.tmp.json)
mod_bam_right_subset=$(jq -r '.mod_bam_right_subset' "$op_dir"/bam_subset_paths.tmp.json)

# remove temporary file
rm "$op_dir"/bam_subset_paths.tmp.json

# create a file that will contain all submitted jobs for plotting (delete if it already exists)
submitted_jobs="$op_dir"/submitted_jobs.txt
rm -rf "$submitted_jobs"

# Define the sbatch_record function that stores the command to be submitted in a file
sbatch_record() {
    # Append the command to submitted_jobs file
    echo "bash $*" >> "$submitted_jobs"
}

# read data line by line from file and make plots
while IFS=' ' read -r -a line_data; do

  read_id="${line_data[0]}"
  feature_start="${line_data[1]:-0}"
  feature_end="${line_data[2]:-0}"

  echo "Plotting read: $read_id with feature start: $feature_start and feature end: $feature_end"

  # strip trailing and leading whitespaces from read id, feature start and feature end
  read_id="$(echo -e "${read_id}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
  feature_start="$(echo -e "${feature_start}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"
  feature_end="$(echo -e "${feature_end}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')"

  # plot analogue-modification and associated data
  sbatch_record plot_one_read_with_features.sh ${plot_option:-} \
      "$read_id" "$mod_bam_subset" "$feature_start" "$feature_end" "$forksense_dir" "$pause_file" "$fa" "$op_dir" \
      "$window_size" "$threshold";

  # plot forksense left
  if [[ -f "$mod_bam_left_subset" ]]; then
    sbatch_record limited_plot_read.sh "$mod_bam_left_subset" "$read_id" 0 "$op_dir" \
      forkSense_left;
  fi

  # plot forksense right
  if [[ -f "$mod_bam_right_subset" ]]; then
    sbatch_record limited_plot_read.sh "$mod_bam_right_subset" "$read_id" 0 "$op_dir" \
      forkSense_right;
  fi

  # plot query coordinate vs reference coordinate
  if [[ -f "$bam_file_subset" ]]; then
    sbatch_record plot_ref_to_query_from_read.sh "$read_id" "$bam_file_subset" \
      "$op_dir"/forkSense_check_using_ref_query;
  fi

done < "$read_id_processed"

# submit plotting jobs as an array job
n_jobs=$(wc -l < "$submitted_jobs")
cd "$main_dir"
launchJob=$(sbatch -c 1 --mem=10G --time=1:00:00 -p ei-medium --array=1-"$n_jobs" -J plotReadArray --mail-user= \
  -D "$current_dir" job_array_converter.sh "$submitted_jobs")
jid=${launchJob##* }
cd "$current_dir"

echo "All plotting jobs have been submitted."
echo "Recapping plotted reads in alphabetical order: "
sort "$read_id_processed" | awk '{print $1}'

# remove temporary file
rm "$read_id_processed"

# collate plots into a pdf
launchJob=$(sbatch -p ei-medium --time=1:30:00 --mail-user= --dependency=afterany:"$jid" \
  convert_plots_to_latex.sh "$pause_file" "$op_dir" "$forksense_dir" "$op_dir" forkSense features)
jid_latex=${launchJob##* }

# launch a dummy job so that we wait till plot jobs have finished
sbatch --dependency=afterany:"$jid":"$jid_latex" --time=00:02:00 --job-name=sleep_10 \
  --wait -p ei-short --wrap="sleep 10;"

# print status of plotting jobs
echo "Plotting job status"
sacct -j "$jid" --format=JobID,State,ExitCode,MaxRSS,Elapsed,Start,End
echo "Pdf compilation job status"
sacct -j "$jid_latex" --format=JobID,State,ExitCode,MaxRSS,Elapsed,Start,End