#!/bin/bash

# goal
# ====
# Given bed file with some reads that map to a contig that consists solely of repeats,
# convert all coordinates to coordinates of the first repeat unit.

# usage
# =====
# bash convert_normal_bed_coords_to_periodic_bed_coords_rDNA.sh <bed_file> <ref_fasta_index> <contig_name>\
#      <num_repeat_units> <intensive_or_extensive> > <output_bed_file>
# NOTE: \ is used to indicate that the command continues on the next line. No spaces are allowed after the \.
# NOTE: > is used to indicate that the output is redirected to a file
# NOTE: all texts within <> need to be specified by the user
# bed_file: bed file with some reads that map to a contig that consists solely of repeats.
#           must have 6 columns: chrom, chromStart, chromEnd, name, score, strand_and_everything_else.
#           only coordinates of entries that correspond to the contig_name will be converted.
#           all other entries will be ignored.
# ref_fasta_index: fasta index file of the reference genome. must contain contig_name.
#            By index, I mean a .fai file that accompanies the fasta file.
#            If not available, use the following command to create the .fai file: samtools faidx <ref_fasta>
# contig_name: name of the contig that consists solely of repeats.
# num_repeat_units: number of repeat units in the contig of interest.
# output_bed_file: output will be redirected to this file if specified.
#                  The file will be created if it does not exist, and overwritten if it does exist.
#                  otherwise, output will be printed to the terminal.
#                  output is the same as the input bed file, except that the coordinates of entries that
#                  correspond to the contig_name will be converted.
# intensive_or_extensive: is the bed file score column intensive or extensive?
#                         if it is intensive, then any interval that overlaps with the periodic boundary will be split,
#                         and the score of the split intervals will be the same as the original interval.
#                         if it is extensive, then the scores of the split intervals will add up to the original
#                         interval, and will be proportional to the length of each interval.

# example
# =======
# bed file contents:
# chr1    1000    2000    read1   7       +
# chr1    2000    3000    read2   2       +
# chr2    1000    3000    read3   3       -
# chr1    8000    9000    read4   2       -
# let's say that chr1 is the contig of interest, it's of length 10000 and that it consists of 4 repeat units.
# so, the first repeat unit is chr1:0-2500, the second repeat unit is chr1:2500-5000, etc.
# so, the output bed file with the intensive option will look like this (order of entries may be different):
# chr1    0       500     read2   2       +
# chr1    1000    2000    read1   7       +
# chr1    2000    2500    read2   2       +
# chr2    1000    3000    read3   3       -
# chr1    500     1500    read4   2       -
# and the output bed file with the extensive option will look like this (order of entries may be different):
# chr1    0       500     read2   1       +
# chr1    1000    2000    read1   7       +
# chr1    2000    2500    read2   1       +
# chr2    1000    3000    read3   3       -
# chr1    500     1500    read4   2       -
# in both cases, reads 1,2 and 4 were considered and reads 2 and 4 were modified.
# read2 which spanned the periodic boundary was split into two intervals.
# in the intensive case, the score of the split intervals is the same as the original interval.
# in the extensive case, the score of the split intervals is proportional to the length of each interval.

# fail if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load user configuration information
source config.sh

# load python
source load_package.sh -python

# check that the correct number of arguments were provided
if [ $# -ne 5 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash convert_normal_bed_coords_to_periodic_bed_coords_rDNA.sh <bed_file> <ref_fasta_index>\ "
    >&2 echo "            <contig_name> <num_repeat_units> <intensive_or_extensive> > <output_bed_file>"
    >&2 echo "NOTE: \ is used to indicate that the command continues on the next line. No spaces are allowed after the \ "
    >&2 echo "NOTE: > is used to indicate that the output is redirected to a file"
    >&2 echo "NOTE: all texts within <> need to be specified by the user"
    >&2 echo "bed_file: bed file with some reads that map to a contig that consists solely of repeats."
    >&2 echo "          must have 6 columns: chrom, chromStart, chromEnd, name, score, strand_and_everything_else."
    >&2 echo "ref_fasta_index: fasta index file of the reference genome. must contain contig_name."
    >&2 echo "contig_name: name of the contig that consists solely of repeats."
    >&2 echo "num_repeat_units: number of repeat units in the contig of interest."
    >&2 echo "intensive_or_extensive: set to intensive/extensive to indicate whether the bed file score column "
    >&2 echo "                        corresponds to an intensive quantity (like density etc.) or an extensive "
    >&2 echo "                        quantity (like count etc.)"
    >&2 echo "output_bed_file: output will be redirected to this file if specified."
    exit 1;
fi

# assign arguments to variables
bed_file=$1
ref_fasta_index=$2
contig_name=$3
num_repeat_units=$4
intensive_or_extensive=$5

# check that the bed file exists
if [ ! -f "$bed_file" ]; then
    >&2 echo "ERROR: bed file does not exist"
    exit 1;
fi

# check that the bed file is in the correct format
if [ ! "$(< "$bed_file" python validate_bed_format.py --allow-float-score)" == "valid" ]; then
    >&2 echo "Error: $bed_file is not in the correct format."
    exit 1;
fi

# check that the bed file coordinates are valid
if [ ! "$(< "$bed_file" python validate_bed_against_fai.py "$ref_fasta_index" )" == "valid"  ]; then
  >&2 echo "Error: $bed_file does not have valid coordinates."
  exit 1;
fi

# check that the fasta index file exists
if [ ! -f "$ref_fasta_index" ]; then
    >&2 echo "ERROR: fasta index file does not exist"
    exit 1;
fi

# check that the contig name is valid
if ! [ "$(grep -c "$contig_name" "$ref_fasta_index")" -eq 1 ]; then
    >&2 echo "ERROR: contig name is not valid"
    exit 1;
fi

# check that the number of repeat units is positive
if [ ! "$num_repeat_units" -gt 0 ]; then
    >&2 echo "ERROR: number of repeat units is not positive"
    exit 1;
fi

# check that the number of repeat units is an integer
if [ "$num_repeat_units" != "${num_repeat_units%.*}" ]; then
    >&2 echo "ERROR: number of repeat units is not an integer"
    exit 1;
fi

# check that the intensive_or_extensive argument is valid
if [ "$intensive_or_extensive" != "intensive" ] && [ "$intensive_or_extensive" != "extensive" ]; then
    >&2 echo "ERROR: intensive_or_extensive argument is not valid"
    exit 1;
fi

# get the length of the contig of interest and the length of each repeat unit
contig_length=$(< "$ref_fasta_index" grep "$contig_name" | awk '{print $2}')
repeat_length=$((contig_length/num_repeat_units))
repeat_length_remainder=$((contig_length % num_repeat_units))

# check that the contig length is a multiple of the number of repeat units
if [ "$repeat_length_remainder" -ne 0 ]; then
    >&2 echo "ERROR: contig length is not a multiple of the number of repeat units"
    exit 1;
fi

# check that the bed files are in the correct format
if [ ! "$(< "$bed_file" python validate_bed_format.py --allow-float-score --six-columns)" == "valid" ]; then
    >&2 echo "Error: $bed_file is not in the correct format."
    exit 1;
fi

# print the script information
echo "# from commit ${COMMITSTR:-NA} generated at ${TIMENOW:-NA} by ${config[name]:-NA} <${config[email]:-NA}>";
echo "# script: $0";
echo "# arguments: $*";

# loop over all entries in the bed file
while read -r line; do

    # if the current entry is a comment, then print it as is
    if [[ "$line" =~ ^#.* ]]; then
      echo "$line"
      continue
    fi

    # extract the coordinates of the current entry
    chrom=$(echo "$line" | awk '{print $1}')
    start=$(echo "$line" | awk '{print $2}')
    end=$(echo "$line" | awk '{print $3}')
    name=$(echo "$line" | awk '{print $4}')
    score=$(echo "$line" | awk '{print $5}')
    strand_and_everything_else=$(echo "$line" | cut -f 6- )

    # if the current entry is not on the contig of interest, then print it as is
    if [ "$chrom" != "$contig_name" ]; then
        echo "$line"
        continue
    fi

    # check that start and end of the current entry are positive, or print an error and exit
    if [ "$start" -lt 0 ] || [ "$end" -lt 0 ]; then
        >&2 echo "ERROR: start and end of the current entry are not positive"
        exit 1;
    fi

    # if the current entry is on the contig of interest, then convert its coordinates and print it

    # decrement start and end by repeat_length
    # until start is less than repeat_length
    while [ "$start" -ge "$repeat_length" ]; do
        start=$((start-repeat_length))
        end=$((end-repeat_length))
    done

    # if the current entry is entirely within one repeat unit, then print it as is
    if [ "$start" -lt "$repeat_length" ] && \
       [ "$end" -le "$repeat_length" ]; then
        echo -e "$chrom\t$start\t$end\t$name\t$score\t$strand_and_everything_else"
        continue
    fi

    # get the length of the current entry
    read_length=$((end-start))

    # if end is greater than twice repeat_length,
    # then decrement it by repeat_length until it is less than twice repeat_length,
    # while printing the interval corresponding to each length removed
    while [ "$end" -gt "$((2*repeat_length))" ]; do
        # calculate the score of the interval depending on whether the bed file score column
        # corresponds to an intensive or extensive quantity
        if [ "$intensive_or_extensive" == "intensive" ]; then
            interval_score=$score
        elif [ "$intensive_or_extensive" == "extensive" ]; then
            interval_score=$(echo "$score"/"$read_length"*"$repeat_length" | bc -l)
        else
            >&2 echo "ERROR: intensive_or_extensive argument is not valid"
            exit 1;
        fi
        echo -e "$chrom\t0\t$repeat_length\t$name\t$interval_score\t$strand_and_everything_else"
        end=$((end-repeat_length))
    done

    # calculate the score of the two pieces to be output depending on whether the bed file score
    # column corresponds to an intensive or extensive quantity
    if [ "$intensive_or_extensive" == "intensive" ]; then
        interval_score_1=$score
        interval_score_2=$score
    elif [ "$intensive_or_extensive" == "extensive" ]; then
        interval_score_1=$(echo "$score"/"$read_length"*"$((repeat_length-start))" | bc -l)
        interval_score_2=$(echo "$score"/"$read_length"*"$((end-repeat_length))" | bc -l)
    else
        >&2 echo "ERROR: intensive_or_extensive argument is not valid"
        exit 1;
    fi

    # print two entries corresponding to the two pieces the current entry is cut into
    echo -e "$chrom\t$start\t$repeat_length\t$name\t$interval_score_1\t$strand_and_everything_else"
    echo -e "$chrom\t0\t$((end-repeat_length))\t$name\t$interval_score_2\t$strand_and_everything_else"

done < "$bed_file"
