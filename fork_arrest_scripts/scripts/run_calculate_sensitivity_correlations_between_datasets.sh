#!/bin/bash

# goal
# -----
# Given sigmoidal strain datasets, calculate pairwise correlations between the pause sensitivities

# details
# -------
# The background is that we have strains where BrdU decreases versus time when synchronized cells are exposed to a
# constant BrdU concentration. We want to determine if our pause calling ability is similar between strains.
# In our pause detection pipeline (this may not be a part of the pipeline script), we calculate sensitivities
# to pause detection in bins across the genome (most likely 1 kb bins).
# We want to determine if the sensitivities are similar between strains.
# We will do this by calculating pairwise correlations between the sensitivities.

# usage
#------
# bash run_calculate_sensitivity_correlations_between_datasets.sh json_file
# can use sbatch but if bedgraph sizes are reasonable, then a bash job on a single node should be sufficient
# json_file: path to json file whose format is described below

# json file format
# ----------------
# A sample file with three datasets is shown below.
# - You can add as many datasets as you want.
# - At least two datasets are required.
# - If/when making correlation plots, the names of the datasets will be used as the labels,
#   so it is best to keep them short.
# - The names can be strain names, dates, or anything else that is useful to you.
# - Sensitivity files are bedgraphs of the sensitivities to pause detection in bins across the genome.
# - All sensitivity files must have the same genomic bins.
# [
#   {
#     "name": "data_1",
#     "sensitivity": "path to sensitivity file 1"
#   },
#   {
#     "name": "data_2",
#     "sensitivity": "path to sensitivity file 2"
#   },
#   {
#     "name": "data_3",
#     "sensitivity": "path to sensitivity file 3"
#   }
# ]

# output
# ------
# The output is a (tab-separated) matrix of correlations between the datasets written to stdout.
# It looks like the following:
# # some comments
# # some more comments
# name data_1 data_2 data_3
# data_1 1.000000 0.100000 0.200000
# data_2 0.100000 1.000000 0.30000
# data_3 0.200000 0.300000 1.000000
# NOTE: as you see above, the matrix is symmetric, so the upper and lower triangles are the same.
#       The diagonal is all 1's because the correlation of a dataset with itself is 1.
# NOTE: the program does not literally write # some comments and # some more comments.
#       These are just placeholders for actual comments.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -jq -python -miller

# load configuration
source config.sh

# function to set calling script information
insert_calling_script_header() {
  sed '1i'\
'# from commit '"${COMMITSTR:-NA}"' generated at '"${TIMENOW:-NA}"' by '"${config[name]:-NA}"' <'"${config[email]:-NA}"'>\n'\
'# script: '"$0"'\n'\
'# arguments: '"$*"'\n'\
"# slurm job name: ${SLURM_JOB_NAME:-NA}"
}

# assign arguments to variables
json_file=${1:-}

# check that the correct number of arguments were provided
if [ "$#" -lt 1 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash $0 json_file"
    >&2 echo "goal: calculate similarities between datasets/strains using calculated pause sensitivities"
    >&2 echo "For more details on what is in the json file, see the comments in the script."
    exit 1;
fi

# check that json file exists
if [ ! -f "$json_file" ]; then
    >&2 echo "ERROR: json file does not exist"
    >&2 echo "json_file: $json_file"
    exit 1;
fi

# using the package jq, create an array of the names of the datasets and an array of the paths to the sensitivity files
mapfile -t names < <(jq -r '.[].name' "$json_file")
mapfile -t sensitivities < <(jq -r '.[].sensitivity' "$json_file")

# check that the name array has at least two elements
if [ "${#names[@]}" -lt 2 ]; then
    >&2 echo "ERROR: at least two datasets are required"
    >&2 echo "number of datasets: ${#names[@]}"
    exit 1;
fi

# check that the number of names and sensitivities are the same
if [ "${#names[@]}" -ne "${#sensitivities[@]}" ]; then
    >&2 echo "ERROR: the number of names and sensitivities are not the same"
    >&2 echo "number of names: ${#names[@]}"
    >&2 echo "number of sensitivities: ${#sensitivities[@]}"
    exit 1;
fi

# run the python script to calculate the correlations between the datasets, then reshape output into a matrix.
# NOTE: we do calculations twice below, ignoring that the matrix is symmetric, because we are lazy.
{
  echo -e "name\tname_2\tvalue"
  for ((i=0;i<${#names[@]};++i)); do
    for ((j=0;j<${#names[@]};++j)); do
      # get the names of the datasets
      name1=${names[$i]}
      name2=${names[$j]}
      # get the paths to the sensitivity files
      sensitivity1=${sensitivities[$i]}
      sensitivity2=${sensitivities[$j]}
      value=$(python correlation_coefficient_bedgraph_values.py "$sensitivity1" "$sensitivity2")
      echo -e "$name1\t$name2\t$value"
    done
  done
} | mlr --tsv reshape -s name_2,value | insert_calling_script_header "$json_file"