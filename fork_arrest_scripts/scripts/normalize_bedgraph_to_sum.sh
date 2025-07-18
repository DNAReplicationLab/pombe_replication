#!/bin/bash

# goal
# -----
# Given a bedgraph, normalize it so that the fourth column sums to a given value

# usage
#------
# cat file.bedgraph | bash normalize_bedgraph_to_sum.sh value
# or
# < file.bedgraph bash normalize_bedgraph_to_sum.sh value

# outputs
# -------
# A four-column bedgraph file is sent to stdout with the fourth column normalized to sum to the given value.

# stop execution if any command fails
set -e

# load configuration
source config.sh

# set temporary directory and make it
tmpDir="${config[scratchDir]:-}"/tmp/temp$(openssl rand -hex 6)
mkdir -p "$tmpDir"

# assign arguments to variables
desired_value=${1:-}

# check that the correct number of arguments were provided
if [ "$#" -lt 1 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: cat file.bedgraph | bash normalize_bedgraph_to_sum.sh desired_value"
    >&2 echo "alter the fourth column of a bedgraph so that it sums to the desired value"
    >&2 echo "desired_value: desired value for the sum of the fourth column"
    exit 1;
fi

# check that the desired value is a number
if ! [[ "$desired_value" =~ ^[0-9]+\.?[0-9]*$ ]]; then
    >&2 echo "ERROR: desired value must be a number"
    exit 1;
fi

# check that desired value is positive
if [ "$(echo "$desired_value <= 0" | bc -l)" -eq 1 ]; then
    >&2 echo "ERROR: desired value must be positive"
    exit 1;
fi

# receive bedgraph from stdin
tmp_bedgraph_stdin=$(mktemp -p "$tmpDir" tmp_bed_file_stdin_XXXXXX.bed)
cat > "$tmp_bedgraph_stdin"

# first, find the sum of the fourth column of the bedgraph
bedgraph_sum=$(< "$tmp_bedgraph_stdin" grep -E -v '^browser|^track|^#' | awk 'BEGIN{sum=0}{sum+=$4}END{print sum}')

# ensure that bedgraph sum is positive
if [ "$(echo "$bedgraph_sum <= 0" | bc -l)" -eq 1 ]; then
    >&2 echo "ERROR: bedgraph sum must be positive. there's probably something wrong with your bedgraph."
    exit 1;
fi

# output normalized bedgraph to stdout
{
  echo "# normalized to sum to $desired_value using normalize_bedgraph_to_sum.sh on $(date)";
  < "$tmp_bedgraph_stdin" grep -E '^browser|^track|^#';
  < "$tmp_bedgraph_stdin" grep -E -v '^browser|^track|^#' |\
    awk -v desired_value="$desired_value" -v bedgraph_sum="$bedgraph_sum"\
      'BEGIN{factor=desired_value/bedgraph_sum;OFS=" "}{print $1, $2, $3, $4*factor}'
}

# delete temporary files
rm -rf "$tmpDir"