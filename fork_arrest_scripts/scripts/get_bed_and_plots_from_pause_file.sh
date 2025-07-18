#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J getBedAndPlotsFromPauseFile
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# goal
# -----
# Given a pause file, plot individual reads and get bed file of pause sites
# NOTE: the pause file is expected to be of a few lines. Do not run this script on pause files with lots of lines.

# usage
#------
# bash get_bed_and_plots_from_pause_file.sh $pause_file $mod_bam $fasta_file $op_dir $json_options $bam_file
# pause_file: file with pause information. This is a tab separated file with headers.
#             Look at convert_pause_file_to_bed.py for a note about the format.
# mod_bam: the .mod.bam file that contains analogue modification probabilities per read. this file must include the
#          reads of interest in the pause file and must be indexed.
# fasta_file: reference genome. must be indexed.
# op_dir: directory where the plots and bed file will be written, preferably an empty directory.
#         do not use a directory with files in it unless you are sure that nothing will be overwritten.
# json_options: (optional) json file with optional fields. See below for details. Defaults to no file if left blank.
# bam_file: (optional) the bam file with alignment details (CIGAR string etc.) of the reads of interest in the pause
#           file. This file must be indexed. By default, this is not used. If used, we plot query vs ref coordinate,
#           and this can be used to watch out for features like insertions/deletions. If you are setting this to a
#           mod BAM file, then you must check that this is a file with alignment information, because not all mod BAM
#           files have this. For e.g.: mod BAM files generated from .detect files do not have alignment information,
#           as DNAscent detect software removes data that is not on the reference.

# json options
# ------------
# The json options file contains optional fields that can be used to add stuff to the plots.
# One needn't specify the file or any of the fields if one doesn't want to.
# The fields are listed below and the json format is explained in the section "json file format".

# mod_bam_left: path to mod bam file with left fork probabilities
# mod_bam_right: path to mod bam file with right fork probabilities
# forksense_dir: directory with forksense files that contain origin, left fork, right fork, termination information.
# overall_pause_file: file with all the pause information in the dataset.

# json file format
# ----------------
# NOTE: A json file has the format as shown below:
# {
#   "key1": "value1",
#   "key2": value2,
#   ...,
#   "keyN": "valueN"
# }
# NOTE: ... does not literally mean three dots. It means that there can be more keys.
# NOTE: numeric values do not need to be quoted, whereas string values need to be quoted.

# outputs
# -------
# Two bed files and a bunch of plots and associated files are sent to the output folder.
# Bed files are a list of pauses, and the same pauses with flanking regions added
# for easy visualization in a genome browser.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages
source load_package.sh -jq -python -bedtools

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
mod_bam=${2:-}
fasta_file=${3:-}
op_dir=${4:-}
json_options=${5:-/dev/null}
bam_file=${6:-/dev/null}

# check that the correct number of arguments were provided
if [ "$#" -lt 4 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash get_bed_and_plots_from_pause_file.sh \$pause_file \$mod_bam \$fasta_file \$op_dir \$json_options \$bam_file"
    >&2 echo "For more details on what these parameters mean, see the comments in the script."
    >&2 echo "The suffix _options means that the parameter is optional."
    exit 1;
fi

# check that input files exist
if [ ! -f "$pause_file" ]; then
    >&2 echo "ERROR: pause file $pause_file does not exist"
    exit 1;
fi

if [ ! -f "$mod_bam" ]; then
    >&2 echo "ERROR: mod bam file $mod_bam does not exist"
    exit 1;
fi

if [ ! -f "$fasta_file" ]; then
    >&2 echo "ERROR: fasta file $fasta_file does not exist"
    exit 1;
fi

# check that the mod bam file is indexed
if [ ! -f "$mod_bam".bai ]; then
    >&2 echo "ERROR: mod bam file $mod_bam is not indexed"
    exit 1;
fi

# check that the fasta file is indexed
if [ ! -f "$fasta_file".fai ]; then
    >&2 echo "ERROR: fasta file $fasta_file is not indexed"
    exit 1;
fi

# check if bam file and the corresponding index exist
if [ "$bam_file" != "/dev/null" ]; then
    if [ ! -f "$bam_file" ]; then
        >&2 echo "ERROR: bam file $bam_file does not exist"
        exit 1;
    fi

    if [ ! -f "$bam_file".bai ]; then
        >&2 echo "ERROR: bam file $bam_file is not indexed"
        exit 1;
    fi
fi

# make output directory
mkdir -p "$op_dir"
op_dir=$(realpath "$op_dir")

# if json options file is provided, check that it exists, and if it does, set the variables
if [ "$json_options" != "/dev/null" ]; then
    if [ ! -f "$json_options" ]; then
        >&2 echo "ERROR: json options file $json_options does not exist"
        exit 1;
    fi
    # set the variables
    mod_bam_left=$(jq -r '.mod_bam_left' "$json_options")
    mod_bam_right=$(jq -r '.mod_bam_right' "$json_options")
    forksense_dir=$(jq -r '.forksense_dir' "$json_options")
    overall_pause_file=$(jq -r '.overall_pause_file' "$json_options")
else
    # set the variables to null
    mod_bam_left="/dev/null"
    mod_bam_right="/dev/null"
    forksense_dir="/dev/null"

    # set overall pause file to pause_file if not set
    overall_pause_file="$pause_file"
fi

# check that the pause file is valid
if [ ! "$(< "$pause_file" python validate_pause_format.py)" == "valid" ]; then
  >&2 echo "Error: pause file is not valid."
  exit 1;
fi

# if there are more than 100 lines in the pause file, then print an error and exit
if [ "$(wc -l < "$pause_file")" -gt 100 ]; then
  >&2 echo "Error: pause file has more than 100 lines. This is alright but it is recommended to use a smaller file."
  >&2 echo "Like, select a few lines from the pause file to make a smaller file and try again."
  >&2 echo "If you want to continue anyway, then comment out the exit statement in the script."
  exit 1;
fi

# convert pause file to bed in two ways
# - first, do the direct conversion
< "$pause_file" python convert_pause_file_to_bed.py --discardUnkeptPauses  > "$op_dir"/pause_file.bed
# - then, do the direct conversion, but use fork directions for strand and slop by +-500 bp
< "$pause_file" python convert_pause_file_to_bed.py --discardUnkeptPauses --LRtoPlusMinus |\
 bedtools slop -i - -g "$fasta_file".fai -b 500 > "$op_dir"/pause_file_LRtoPlusMinus_slop500.bed

# extract read ids
read_id_file="$op_dir"/read_ids.txt
<  "$op_dir"/pause_file.bed python print_valid_data_bed_lines.py | awk '{print $4}' |\
  sort | uniq > "$read_id_file"

# plot the reads
cd plotting_and_short_analyses;
bash plot_n_reads.sh "$read_id_file" "$mod_bam" "$mod_bam_left" "$mod_bam_right" "$forksense_dir" "$overall_pause_file"\
  "$fasta_file" "$op_dir" "$bam_file"