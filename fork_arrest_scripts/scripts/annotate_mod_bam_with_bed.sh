#!/bin/bash
#SBATCH --mem=20G
#SBATCH -c 16
#SBATCH -p ei-medium
#SBATCH -J annotate_mod_bam_with_bed
#SBATCH --mail-type=END,FAIL
#SBATCH --time=2:59:59
#SBATCH --constraint=""

# This script annotates a BAM file with a BED file. It uses the bedtools tag command.
# It's a thin wrapper around bedtools tag, but it does a few things:
# - It standardizes the tag we use for the annotation ('XT').
# - It collects documentation and instructions about the annotation.

# Usage: sbatch annotate_mod_bam_with_bed.sh <input_bam> <output_bam> <bed_file_1>:<label_1> <bed_file_2>:<label_2> ...
#        where ... means you can have as many bed files as you want.
#        <input_bam>: The input BAM file.
#        <output_bam>: The output BAM file with tags for each bed file you provide.
#        <bed_file_1>: The first bed file to use for annotation.
#        <label_1>: The label to use for the tag for the first bed file.

# Example
# =======
# sbatch annotate_mod_bam_with_bed.sh input.bam output.bam sample_a.bed:sample_a sample_b.bed:sample_b

#  If the input bed files are in BED6 format i.e. six tab-separated columns of contig, start, end, name, score,
#  strand, then the output BAM file will have two tags per (intersecting) read that look like
#  XT:Z:sample_a:contig_a1:start_a1-end_a1,name_a1,score_a1,strand_a1,sample_a:contig_a2:\
#  start_a2-end_a2,name_a2,score_a2,strand_a2,...;sample_b:contig_b1:start_b1-end_b1,name_b1,\
#  score_b1,strand_b1,sample_b:contig_b2:start_b2-end_b2,name_b2,score_b2,strand_b2,...
#  NOTE: \ and ... are used to indicate line breaks and "and so on".
#  The actual output will not have a line break or these symbols.

#  If the input bed files are in BED3 format i.e. three tab-separated columns of contig, start, end, then the
#  output BAM file will look similar to the above but the name, score, and strand entries will be empty
#  e.g. XT:Z:sample_a:contig_a1:start_a1-end_a1,,,,sample_a:contig_a2:start_a2-end_a2,,,,...

# fail upon error
set -e

# check that there are at least 3 arguments
if [ $# -lt 3 ]; then
    >&2 echo "Usage: sbatch annotate_mod_bam_with_bed.sh <input_bam> <output_bam> <bed_file_1>:<label_1> <bed_file_2>:<label_2> ..."
    >&2 echo "       where ... means you can have as many bed files as you want."
    >&2 echo "       <input_bam>: The input BAM file."
    >&2 echo "       <output_bam>: The output BAM file with have tags for each bed file you provide."
    >&2 echo "       <bed_file_1>: The first bed file to use for annotation."
    >&2 echo "       <label_1>: The label to use for the tag for the first bed file."
    >&2 echo "       and so on..."
    exit 1;
fi

# check that the first argument is a file that exists
if [ ! -f "$1" ]; then
    >&2 echo "Error: $1 is not a file."
    exit 1;
fi

# collect argument from the third one to the end and store them in an array
bed_files=()
labels=()
for arg in "${@:3}"; do
    bed_file=$(echo "$arg" | cut -d':' -f1)
    label=$(echo "$arg" | cut -d':' -f2)
    bed_files+=("$bed_file")
    labels+=("$label")
done

# check that the bed files exist
for bed_file in "${bed_files[@]}"; do
    if [ ! -f "$bed_file" ]; then
        >&2 echo "Error: $bed_file is not a file."
        exit 1;
    fi
done

# check that none of the labels are empty
for label in "${labels[@]}"; do
    if [ -z "$label" ]; then
        >&2 echo "Error: one of the labels is empty."
        exit 1;
    fi
done

# check that none of the labels have spaces in them
for label in "${labels[@]}"; do
    if [[ "$label" == *" "* ]]; then
        >&2 echo "Error: one of the labels has a space in it."
        exit 1;
    fi
done

# check that none of the bed files have spaces in them
for bed_file in "${bed_files[@]}"; do
    if [[ "$bed_file" == *" "* ]]; then
        >&2 echo "Error: one of the bed files has a space in it."
        exit 1;
    fi
done

# create a string of bed files, with a space between each
bed_files_string=$(printf " %s" "${bed_files[@]}")

# create a string of labels, with a space between each
labels_string=$(printf " %s" "${labels[@]}")

# load bedtools
source load_package.sh -bedtools -python -samtools

# load git repo labels
source load_git_repo_labels.sh

# load user configuration variables
source config.sh

# check that the bed files are in the correct format
for bed_file in "${bed_files[@]}"; do
    if [ ! "$(< "$bed_file" python validate_bed_format.py --allow-float-score)" == "valid" ]; then
        >&2 echo "Error: $bed_file is not in the correct format."
        exit 1;
    fi
done

# make a temporary file to store bed data
tmp_file=$(mktemp)

# create the output
# shellcheck disable=SC2086
bedtools tag -i "$1" -files $bed_files_string -labels $labels_string -intervals -tag XT > "$tmp_file"

# create a header for the output
script_name_and_time=$(basename "$0")_${TIMENOW:-NA}
comment_str_1="# $script_name_and_time - from commit ${COMMITSTR:-NA} generated at ${TIMENOW:-NA} by ${config[name]:-NA} <${config[email]:-NA}>"
comment_str_2="# $script_name_and_time - script: $0";
comment_str_3="# $script_name_and_time - arguments: $*";
program_str="bedtools tag -i ""$1"" -files $bed_files_string -labels $labels_string -intervals -tag XT > ""$tmp_file"
bedtools_version=$(bedtools --version | awk '{print $2}')
samtools view -H "$tmp_file" |\
 sed '1i@CO\t'"$comment_str_1" |\
 sed '1i@CO\t'"$comment_str_2" |\
 sed '1i@CO\t'"$comment_str_3" |\
 python add_program_name_to_sam_header.py bedtools "$bedtools_version" "$program_str" > "$tmp_file".header

# create the new header
samtools reheader "$tmp_file".header "$tmp_file" > "$tmp_file".tmp

# sort and index the bam file
samtools sort -@ 16 -o "$2" "$tmp_file".tmp
samtools index "$2"

# remove the temporary files
rm "$tmp_file" "$tmp_file".header "$tmp_file".tmp