#!/bin/bash

# goal
# -----
# After obtaining pause enrichment across features, compare this result across datasets/nanopore runs.

# usage
#------
# bash run_collect_plots_per_dataset_and_feature_to_latex.sh input_file output_dir
# input_file: tab-separated with columns: dataset, feature, folder
#             First non-comment line must have column names and comments must start with #
#             Folder contains png files.
# output_dir: directory where the output pdf will be saved

# input file format explanation
# -----------------------------
#  e.g. (imagine spaces are tabs):
# dataset feature folder
# johndoe_dataset1 tRNA /path_1/to/pause_enrichment_plots
# janedoe_dataset1 tRNA /path_2/to/pause_enrichment_plots

# in scenario above, john and jane doe have done nanopore runs and we've plotted pauses across tRNAs
# and stored it in various png files in the two folder.
# now, we want to compare the pause enrichment between the two datasets.

# outputs
# -------
# Pdf file to the output folder.
# In the example above, the pdf file will have a section called tRNA, within which subsections will be names of
# png files and each subsubsection within that subsection will be a dataset name and the corresponding png file.
# The goal is you go to a subsection, and compare the subsubsections to see if there are any differences between
# the datasets.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -python -latex

# set temporary directory
tmpDir="${config[scratchDir]:-}"/tmp/temp$(openssl rand -hex 6)
mkdir -p "$tmpDir"


# assign arguments to variables
input_file=${1:-}
output_dir=${2:-}

# check that the correct number of arguments were provided
if [ "$#" -lt 2 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash run_collect_plots_per_dataset_and_feature_to_latex.sh input_file output_dir"
    >&2 echo "input_file: tab-separated with columns: dataset, feature, folder"
    >&2 echo "            First non-comment line must have column names and comments must start with #"
    >&2 echo "            For details, see the comments in this script"
    >&2 echo "output_dir: directory where the output pdf will be saved"
    exit 1;
fi

# check that the input file exists
if [ ! -f "$input_file" ]; then
    >&2 echo "ERROR: input file does not exist"
    >&2 echo "input_file: $input_file"
    exit 1;
fi

# make the output directory if it does not exist
mkdir -p "$output_dir"
output_dir=$(realpath "$output_dir")

# get a temporary file name
sample_tex=$(mktemp --tmpdir="$tmpDir" sample.XXXXXX.tex)

# run the script and make a pdf
common_prefix=collect_pause_plots
{
  grep '^#' "$input_file" | cat
  < "$input_file" python get_png_in_directories.py
}|\
   tee "$output_dir"/"$common_prefix".tsv |\
   python collect_plots_per_dataset_and_feature_to_latex.py > "$sample_tex"
pdflatex -output-directory="$output_dir" -jobname="$common_prefix" "$sample_tex"
pdflatex -output-directory="$output_dir" -jobname="$common_prefix" "$sample_tex"

# clean up
rm -rf "$tmpDir"