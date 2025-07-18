#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J plotPauseProfileFromReportSplitBed
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# goal
# -----
# Background- we have a directory where we split a bed file by length along entries and by some other value
# e.g. transcription rates. In this directory, we have a bunch of json files *pauses_report.json. We want
# to gather the pause statistics from these files and plot them to show pause count vs length,
# pause count vs the value (which could be transcription rate) etc.

# usage
#------
# sbatch run_plot_pause_profiles_from_pause_reports_from_split_bed.sh $pause_report_directory $plot_output_directory
# can use bash instead of sbatch if need be, think about running times and decide accordingly.
# $pause_report_directory: the directory where the *pauses_report.json files are located
# $plot_output_directory: the directory where the plots will be outputted.
#                         Preferably, this should be a new and/or empty directory.
#                         Directory will be created if it doesn't exist.

# outputs
# -------
# A bunch of plots are sent to the specified output directory.
# We write some messages to the standard output/SBATCH output for logging purposes.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -python -miller -jq -R

# load configuration
source config.sh

# assign arguments to variables
pause_report_directory=${1:-}
plot_output_directory=${2:-}

# check that the correct number of arguments were provided
if [ "$#" -lt 2 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: sbatch plot_pause_profiles_from_pause_reports_from_split_bed.sh \$pause_report_directory \$plot_output_directory"
    >&2 echo "\$pause_report_directory: the directory where the *pauses_report.json files are located"
    >&2 echo "\$plot_output_directory: the directory where the plots will be outputted."
    exit 1;
fi

# check that the pause report directory exists
if [ ! -d "$pause_report_directory" ]; then
    >&2 echo "ERROR: pause report directory does not exist"
    >&2 echo "pause report directory: $pause_report_directory"
    exit 1;
fi

# check that the pause report directory contains at least one file called *pauses_report.json
FILES_EXIST=false
for file in "$pause_report_directory"/*pauses_report.json
do
  if [ -e "$file" ]
  then
    FILES_EXIST=true
    break
  fi
done

if [ "$FILES_EXIST" = false ]; then
    >&2 echo "ERROR: pause report directory does not contain any *pauses_report.json files"
    >&2 echo "pause report directory: $pause_report_directory"
    exit 1;
fi

# make the plot output directory if it doesn't exist
mkdir -p "$plot_output_directory"

# make a temporary directory
tmp_dir="$plot_output_directory"/temp_$(openssl rand -hex 6)_made_on_$(date +%Y%m%d)_delete_later
mkdir -p "$tmp_dir"

# make a summary file
echo "{\"directories\": [\"$pause_report_directory\"]}" |\
  bash collate_pause_report_statistics.sh stdin  |\
    mlr -s reorganize_collated_pause_report_statistics.mlr --ofmt '%.1f' --tsv \
      > "$pause_report_directory"/pause_summary.tsv

# validate summary file
if [ ! "$(< "$pause_report_directory"/pause_summary.tsv python validate_summary_file.py)" == "valid" ]; then
    >&2 echo "ERROR: pause summary file is not valid"
    >&2 echo "pause summary file: $pause_report_directory/pause_summary.tsv"
    exit 1;
fi

# log that the summary file was prepared and validated
date +\[%Y_%m_%d_%H:%M:%S\]
echo "Prepared and validated a summary file $pause_report_directory/pause_summary.tsv"

# create a list of relative directions
relative_directions=("all" "co-directional" "head-on")

# log message
date +\[%Y_%m_%d_%H:%M:%S\]
echo "Plotting pause count summed across length for all, co-directional, and head-on pauses, whichever are available"

# We set up plot options for each subset, assuming there are 100 value splits.
# In reality, there may be fewer splits, but we can handle that.
for ((i=-1;i<100;++i));
  do
  if [ "$i" -eq -1 ]; then
    query_str="NA"
    suffix_str=""
  else
    query_str="VS_${i}"
    suffix_str="_VS_${i}"
  fi
  temp_file=$(mktemp --tmpdir="$tmp_dir" temp.XXXXXX.tsv)
  # shellcheck disable=SC2016
  # shellcheck disable=SC1010
  < "$pause_report_directory"/pause_summary.tsv \
    mlr --tsv filter '$division == "'"$query_str"'"' then put '$division = $relative_direction' > "$temp_file"
  if [ "$(wc -l < "$temp_file")" -ne 0 ]; then
    echo "Plotting pause count summed across length for $query_str for all, co-directional, and head-on pauses, whichever are available"
    file_name_prefix="$plot_output_directory"/pause_summary"$suffix_str"
    < "$temp_file" tee "$file_name_prefix".data |\
      Rscript plotting_and_short_analyses/plot_pause_count_across_ranked_genomic_features.R \
        "$file_name_prefix".png \
        "Division" "" "" "" "" "" "" 0.7
  fi
done

# plot all, co-directional and head-on separately for different subsets, setting plot options appropriately
subsets=("^LS_[0-9]+$" "^VS_[0-9]+$")
suffixes=("" "")
prefixes=("length_split" "value_split")
widths=(9 5)
x_label_scales=(0.5 "")
error_bar_whisker_scales=(0.7 "")

# We set up plot options for each subset, assuming there are 100 value splits.
# In reality, there may be fewer splits, but we can handle that.
for ((i=0;i<100;++i));
do
  subsets+=( "^VS_${i}_LS_[0-9]+" )
  suffixes+=( "_VS_${i}" )
  prefixes+=( "length_split" )
  widths+=( 9 )
  x_label_scales+=( 0.5 )
  error_bar_whisker_scales+=( 0.7 )
done

for ((j=0;j<${#subsets[@]};++j));
do

  query_str="${subsets[$j]}"
  suffix_str="${suffixes[$j]}"
  prefix_str="${prefixes[$j]}"

  # - plot bar chart for length split i.e. a meta analysis where we are
  #   looking for pause profile along length of that entry
  # - do this again for co-directional and head-on separately if these are available

  for ((i=0;i<${#relative_directions[@]};++i));
  do

    if [ "${relative_directions[$i]}" == "all" ]; then
      relative_direction_suffix="all"
    elif [ "${relative_directions[$i]}" == "co-directional" ]; then
      relative_direction_suffix="cd"
    elif [ "${relative_directions[$i]}" == "head-on" ]; then
      relative_direction_suffix="ho"
    else
      >&2 echo "ERROR: relative direction not recognized"
      >&2 echo "relative direction: ${relative_directions[$i]}"
      exit 1;
    fi

    temp_file=$(mktemp --tmpdir="$tmp_dir" temp.XXXXXX.tsv)
    # shellcheck disable=SC2016
    # shellcheck disable=SC1010
    < "$pause_report_directory"/pause_summary.tsv \
    mlr --tsv filter '$division =~ "^'"$query_str"'" && $relative_direction=="'"${relative_directions[$i]}"'"'\
      then sort -t division > "$temp_file"

    # if temp file is empty, skip plotting
    if [ "$(wc -l < "$temp_file")" -eq 0 ]; then
      continue
    fi

    # log message
    date +\[%Y_%m_%d_%H:%M:%S\]
    echo "Plotting pause count profile for $query_str and ${relative_directions[$i]}"

    file_name_prefix="$plot_output_directory"/"$prefix_str"_"${relative_direction_suffix}""$suffix_str"
    < "$temp_file" tee "$file_name_prefix".data |\
      Rscript plotting_and_short_analyses/plot_pause_count_across_ranked_genomic_features.R\
        "$file_name_prefix".png "Division" ""\
          "${widths[$j]}" "" "${x_label_scales[$j]}" "" "" 0.7 "${error_bar_whisker_scales[$j]}"

  done

done

rm -rf "$tmp_dir"