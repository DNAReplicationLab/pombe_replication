#!/bin/bash

# goal
# -----
# Given a bed file with numerical values in a column, split bed file based on descending value and segments along each
# bed entry (see example below).

# example
# -------
# let's say the input bed file is (replace spaces with tabs):
# chr1  100  200  A  1  + 7
# chr2  400  500  B  10  -  5
# chr3  500  600  A  1  + 12
# chr4  800  900  B  10  -  15
# and the user requests 2 splits according to the seventh column,
# split by 50 bp segments in a head-to-tail manner,
# and requests the prefix file to be "file".
# Then eight files are produced, with the ninth just being a repeat of the input file.
# (VS means value split, LS means length split).
# file_VS_1.bed:
# chr3  500  600  A  1  + 12
# chr4  800  900  B  10  -  15
# file_VS_1_LS_1.bed:
# chr3  550  600  A  1  + 12
# chr4  800  850  B  10  -  15
# file_VS_1_LS_2.bed:
# chr3  500  550  A  1  + 12
# chr4  850  900  B  10  -  15
# file_VS_2.bed:
# chr1  100  200  A  1  + 7
# chr2  400  500  B  10  -  5
# file_VS_2_LS_1.bed:
# chr1  150  200  A  1  + 7
# chr2  400  450  B  10  -  5
# file_VS_2_LS_2.bed:
# chr1  100  150  A  1  + 7
# chr2  450  500  B  10  -  5
# input bed file:
# chr1  100  200  A  1  + 7
# chr2  400  500  B  10  -  5
# chr3  500  600  A  1  + 12
# chr4  800  900  B  10  -  15
# file_LS_1.bed:
# chr1  150  200  A  1  + 7
# chr2  400  450  B  10  -  5
# chr3  550  600  A  1  + 12
# chr4  800  850  B  10  -  15
# file_LS_2.bed:
# chr1  100  150  A  1  + 7
# chr2  450  500  B  10  -  5
# chr3  500  550  A  1  + 12
# chr4  850  900  B  10  -  15

# in algorithm v2, similar inputs and outputs are required, but with a few changes:
# - if you request split by length, the input bed file must have zero interval sizes, and we grow the intervals
#   ourselves by +- so many bp according to another user-requested parameter.
# - w.r.t output, overlaps between the length-split files are removed. So any position on the genome will only
#   appear in one file among the series LS_1, LS_2, ... etc. This is not the case in v1. This is not evident
#   in the example above as the intervals are all far apart.

# usage
#------
# bash split_bed_file_by_column_value_and_along_length.sh input.json
# where input.json is a json file with the following parameters (note how numbers or boolean values are not within
# quotes):
# {
#   "input_bed_file": "/path/to/bed/file.bed",
#   "column_of_interest": 7,
#   "num_split_by_value": 2,
#   "do_equal_split_by_genome_coverage_not_number_of_entries": false,
#   "split_by_length_bp": 50,
#   "align_by": "head-to-tail",
#   "output_prefix": "/path/to/output/directory/prefix_of_file",
#   "feature": "tRNA",
#   "split_algorithm_version": "v1",
#   "grow_region_bp_if_version_is_v2": 100,
#   "fai_file": "/path/to/fasta/file.fasta.fai"
# }
# Meanings of parameters:
# input_bed_file: input bed file
# column_of_interest: (default -1 i.e. invalid/unused) (1-based) column of interest to sort in descending order and
#                     split by. must be numerical column. 1-based means 1 is the first column.
#                     If you do not set this but set num_split_by_value > 1, you will get an error.
# num_split_by_value: resultant number of split files using the column of interest. If you want to split into 2 files,
#                     use 2. If you want no splits, use 1.
# do_equal_split_by_genome_coverage_not_number_of_entries: (default false) if true, then bed file is split such that
#                                                          each piece has (approx) equal genome coverage.
#                                                          Default is false, which means that each piece has (approx)
#                                                          equal number of entries.
#                                                          What you use depends on your application.
#                                                          If there's no correlation between bed entry length and the
#                                                          value of interest, then it should make almost no difference
#                                                          what you use here.
# split_by_length_bp: number of base pairs to split each bed entry by. If you want no splits, use 0.
# align_by: "head-to-tail" or "tail-to-head" or "increasing-ref" or "decreasing-ref". Default is "head-to-tail".
#           Should be ignored by script (either this one or downstream tools) if split_by_length_bp is 0 as if you are
#           not splitting by length, then there's no need to align.
# output_prefix: prefix to use for output files
# feature: feature in the input bed file e.g. "tRNA", "gene"
# split_algorithm_version: (default "v1") version of the split algorithm to use.
#                          "v1" is the original algorithm and is chosen by default.
#                          "v2" is a new algorithm that grows regions by +- so many bp if requested,
#                           names files differently, and trims some kinds of overlaps.
#                           For further details, see the section "Versions of splitting algorithm" below.
# grow_region_bp_if_version_is_v2: (default 0) regions are grown by +- bp if requested. Only used if version is "v2".
# fai_file: (default "" i.e. unused, required if v2) fasta index file of the reference genome.
#           Only used if version is "v2".

# Versions of splitting algorithm
# -------------------------------
# There are two versions of the splitting algorithm: v1 and v2.
# ==================
# Major difference 1
# ==================
# The major difference in output boils down to how overlaps are handled.
# Let us consider an example: If we have a feature chr1:1000-3000 and another feature chr1:2000-4000 on the
# same + strand and request a tile size of 100 bp,the genomic tile chr1:2200-2300 will appear in bin 13 due to the first
# feature and bin 3 due to the second feature (I am numbering bins from 1 in this example).
# In v1, this is tolerated, whereas in v2, the tile appears only in bin 3 as 3 is lower than 13, so we also
# take into account the closeness of a tile to a feature.
# This difference is irrelevant if there are no overlaps e.g. if the bed file is of centromere positions.
# ==================
# Major difference 2
# ==================
# We grow regions by +- so many bp according to user request, whereas in v1 it is expected that you grow the
# regions yourself before running the script. Accordingly, in v2, if a length split is requested by the user
# (by setting split_by_length_bp in the input json to a positive integer), we only accept bed files where all interval
# sizes are 0. This gives us greater control over the output as we have an even interval size and we can deal
# with edge cases like features at the start or end of a contig more easily as we do the growing ourselves.
# We do not expect the edge cases to affect the output much, but it is a nice feature to have.
# We also expect the grow_region_bp_if_version_is_v2 parameter to be a multiple of the split_by_length_bp
# so that we get clean splits.

# outputs
# -------
# Output is to stdout and is a json object with the list of files created and an identifier for each.
# A list with two sample entries is shown below:
# [
#  {
#    "division": "LS_1",
#    "bed_file": "/blah/file_LS_1.bed",
#    "feature": "bg"
#  },
#  {
#    "division": "LS_2",
#    "bed_file": "/blah/file_LS_2.bed",
#    "feature": "bg"
#  }
# ]

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -python -jq

# load configuration
source config.sh

# check that the correct number of arguments were provided
if [ "$#" -ne 1 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash split_bed_file_by_column_value_and_along_length.sh input.json"
    >&2 echo "where input.json is a json file with the following parameters:"
    >&2 echo "input_bed_file, column_of_interest, num_split_by_value, split_by_length_bp, align_by, output_prefix, feature."
    >&2 echo "For more details, see the comments in the script."
    exit 1;
fi

# read arguments
input_json=${1:-}

# check that the input json file exists
if [ ! -f "$input_json" ]; then
    >&2 echo "ERROR: input json file $input_json does not exist"
    exit 1;
fi

# assign arguments to variables
bed_file=$(jq -r '.input_bed_file' "$input_json")
value_column=$(jq -r '.column_of_interest // "-1"' "$input_json") # default to -1 if not provided i.e. invalid
num_value_split=$(jq -r '.num_split_by_value' "$input_json")
is_equal_split_by_genome_coverage=$(jq -r '.do_equal_split_by_genome_coverage_not_number_of_entries' "$input_json")
length_split=$(jq -r '.split_by_length_bp' "$input_json")
align_by=$(jq -r '.align_by // "head-to-tail"' "$input_json")
op_prefix=$(jq -r '.output_prefix' "$input_json")
feature=$(jq -r '.feature' "$input_json")
split_algorithm_version=$(jq -r '.split_algorithm_version // "v1"' "$input_json")
grow_region_bp_if_version_is_v2=$(jq -r '.grow_region_bp_if_version_is_v2 // 0' "$input_json")
fai_file=$(jq -r '.fai_file // ""' "$input_json")

# ensure that the bed_file exists
if [ ! -f "$bed_file" ]; then
    >&2 echo "ERROR: bed file $bed_file does not exist"
    exit 1;
fi

# ensure that the bed file is valid
if [ ! "$(< "$bed_file" python validate_bed_format.py --six-columns --allow-float-score --no-dot-strand)" == "valid" ];
then
  >&2 echo "Error: bed file needs at least six columns and must not have dot strand."
  exit 1;
fi

# ensure that length_split is a positive integer
if [ -z "$length_split" ] || [[ ! "$length_split" =~ ^[0-9]+$ ]]; then
    >&2 echo "ERROR: length_split must be numeric"
    exit 1;
fi

if [ "$length_split" -lt 0 ]; then
    >&2 echo "ERROR: length_split must be an integer >= 0"
    exit 1;
fi

# ensure that num_value_split is a positive integer
if [ -z "$num_value_split" ] || [[ ! "$num_value_split" =~ ^[0-9]+$ ]]; then
    >&2 echo "ERROR: num_split_by_value must be numeric"
    exit 1;
fi

if [ "$num_value_split" -le 0 ]; then
    >&2 echo "ERROR: num_split_by_value must be an integer > 0"
    exit 1;
fi

# create many bed files as requested, and output a json object with the list of files created and an identifier for each
# about what kind of division was performed called "division"
if [ "$split_algorithm_version" == "v1" ]; then

  {

    # output original bed file
    original_bed_file="[\n{\"bed_file\": \"$bed_file\", \"division\": \"NA\"}\n]"
    echo -e "$original_bed_file"

    # split by value in column requested according to number of pieces requested, if requested
    if [ "$num_value_split" -gt 1 ]; then
      if [ "$is_equal_split_by_genome_coverage" == "true" ]; then
        value_split_bed_files=$(\
          < "$bed_file" python split_bed_file_by_column_value_desc.py --prefix "$op_prefix"_VS\
            --column "$value_column" --num-files "$num_value_split" --equal-length-instead-of-equal-rows |\
              jq 'map(.division |= "VS_\(.)")'\
          )
      else
        value_split_bed_files=$(\
        < "$bed_file" python split_bed_file_by_column_value_desc.py --prefix "$op_prefix"_VS\
          --column "$value_column" --num-files "$num_value_split" | jq 'map(.division |= "VS_\(.)")'\
          )
      fi
      echo "$value_split_bed_files"
    fi

    # split each bed file entry by length requested, if requested
    if [ "$length_split" -gt 0 ]; then
      < "$bed_file" python split_each_bed_file_entry_along_length.py --prefix "$op_prefix"_LS\
        --num-bases "$length_split" --align-by "$align_by" | jq 'map(.division |= "LS_\(.)")'

      # now, take each split-by-value bed file and split by length
      if [ "$num_value_split" -gt 1 ]; then
        division_array=()
        while IFS='' read -r line; do division_array+=("$line"); done < <(echo -e "$value_split_bed_files" | jq -r '.[].division')
        bed_files_array=()
        while IFS='' read -r line; do bed_files_array+=("$line"); done < <(echo -e "$value_split_bed_files" | jq -r '.[].bed_file')

        for i in "${!division_array[@]}"; do
          division_i="${division_array[i]}"
          bed_file_i="${bed_files_array[i]}"
          < "$bed_file_i" python split_each_bed_file_entry_along_length.py --prefix "$op_prefix"_"$division_i"_LS\
            --num-bases "$length_split" --align-by "$align_by" | jq 'map(.division |= "'"$division_i"'_LS_\(.)")'
        done
      fi
    fi

  }  | jq '.[]' | jq -s 'sort_by(.division)' | jq '.[].feature = "'"$feature"'"'

elif [ "$split_algorithm_version" == "v2" ]; then

  # check that the fai file exists
  if [ ! -f "$fai_file" ]; then
    >&2 echo "ERROR: fai file $fai_file does not exist"
    exit 1;
  fi

  # set command option for equal genome coverage
  command_option_equal_genome_coverage=""
  if [ "$is_equal_split_by_genome_coverage" == "true" ]; then
    command_option_equal_genome_coverage="--is-value-split-equal-genome-instead-of-equal-rows"
  fi

  # if grow region bp is set, ensure it is >= 0
  if [ "$grow_region_bp_if_version_is_v2" -lt 0 ]; then
    >&2 echo "ERROR: grow_region_bp_if_version_is_v2 must be an integer >= 0"
    exit 1;
  fi

  # if grow region bp is set and a length split is requested,
  # ensure grow region bp is a multiple of the length split,
  # and ensure that bed files have all interval sizes of 0
  if [ "$grow_region_bp_if_version_is_v2" -gt 0 ] && [ "$length_split" -gt 0 ]; then
    if [ "$((grow_region_bp_if_version_is_v2 % length_split))" -ne 0 ]; then
      >&2 echo "ERROR: grow_region_bp_if_version_is_v2 must be a multiple of the length split"
      exit 1;
    fi
    if [ ! "$(< "$bed_file" python validate_bed_format.py --six-columns --allow-float-score --no-dot-strand --zero-base)" == "valid" ]; then
      >&2 echo "ERROR: bed file must have all interval sizes of 0 if grow_region_bp_if_version_is_v2 is set"
      exit 1;
    fi
  fi

  # run the script to split the bed file
  < "$bed_file" \
    python grow_and_split_bed_file_by_column_value_and_length.py \
      --grow-region-num-bases "$grow_region_bp_if_version_is_v2" \
      --value-column "$value_column" \
      --num-value-splits "$num_value_split" \
      --length-split-num-bases "$length_split" \
      --align-by "$align_by" \
      --output-prefix "$op_prefix" \
      --fai-file "$fai_file" \
      $command_option_equal_genome_coverage  | jq '.[].feature = "'"$feature"'"'

else

  >&2 echo "ERROR: split_algorithm_version must be either v1 or v2"
  exit 1;

fi