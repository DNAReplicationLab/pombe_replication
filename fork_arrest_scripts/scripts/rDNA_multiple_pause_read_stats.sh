#!/bin/bash

# goal
# -----
# Collect reads with multiple pauses on them, plot c.d.f. of distances between pauses.

# further details on goal
# ------------------------
# We want to see if there is a correlation between neighbouring rDNA units.
# NOTE: we know that we do not have uniform coverage over all possible neighbouring distances i.e.
#       our coverage on neighbouring units is much higher than our coverage on ten units away due to our read
#       length distributions.

# usage
#------
# bash <script_name.sh> pause_file output_plot
# pause_file is pauses found in the rDNA in our usual format. It's produced by our rDNA scripts and will have a
#   name like pauseList_processed.
# output_plot is the /path/to/output_plot.png

# outputs
# -------
# Two files: one plot with c.d.f. of distances between pauses, and an associated table with the data.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages
source load_package.sh -python -miller -R

# load configuration
source config.sh

# assign arguments to variables
pause_file=${1:-/dev/null}
output_plot=${2:-/dev/null}

# check that the correct number of arguments were provided
if [ "$#" -lt 2 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash <script_name.sh> pause_file output_plot"
    >&2 echo "pause_file is pauses found in the rDNA in our usual format."
    >&2 echo "output_plot is the /path/to/output_plot.png"
    >&2 echo "For further details, see the script header."
    exit 1;
fi

# make output directory if it does not exist
output_dir=$(realpath "$(dirname "$output_plot")")
mkdir -p "$output_dir"

# convert output_plot to absolute path
output_plot=$(realpath "$output_plot")

# check that the pause file exists and is valid
if [ ! -f "$pause_file" ]; then
  >&2 echo "ERROR: $pause_file is not a file or does not exist."
  exit 1
fi
if [ ! "$(< "$pause_file" python validate_pause_format.py)" == "valid" ]; then
  >&2 echo "Error: pause file is not valid."
  exit 1;
fi

# get the distances between pauses, retaining only rightward pauses
< "$pause_file" python convert_pause_file_to_bed.py --discardUnkeptPauses --LRtoPlusMinus |\
  grep -E -v '^#|^browser|^track' |  awk '$6=="+"{print $4, $2}' |\
  sort -k 1,1 -k2,2n | uniq -w 36 --all-repeated=separate  |\
  awk '
	{
		if ($2 != "" && prev != "") {
			diff = $2 - prev
			print diff
		}
		prev = $2
	}' | sed '1idist' | mlr --tsv histogram -f dist --lo 0 --hi 200000 --nbins 200 > "$output_plot".hist

# form the plot
cd plots_for_publication
< "$output_plot".hist Rscript plot_rDNA_multiple_pause_read_stats.R "$output_plot"