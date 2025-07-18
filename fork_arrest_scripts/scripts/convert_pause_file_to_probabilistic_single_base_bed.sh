#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J pauseProbSglBase
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# goal
# -----
# given a pause file obtained by the cut-and-align procedure using a reference sigmoid,
# - extract probability distribution function (p.d.f.) of distance of pause site from inflection point per sigmoid
# - on non-paused forks, calculate probability of pause detection in fork using distance from inflection point and pdf.

# usage
#------
# sbatch convert_pause_file_to_probabilistic_single_base_bed.sh $pause_file $use_only_these_keep_columns $fasta $decimate_T_optional
# pause_file: pause file obtained by the cut-and-align procedure using a reference sigmoid. This is in tab-separated
#             format with a header line with column names, comments starting with '#' and
#             required columns detectIndex, pauseSite, paramStrNoPause, leftPieceXmax, rightPieceXmin, width and
#             optional columns keep*.
#             - detectIndex is a string of the format readid_contig_alignStart_alignEnd_alignOrn_forkDirn_forkStart_forkEnd.
#               names are self-explanatory. alignOrn = fwd/rev and forkDirn = L/R.
#               alignStart < alignEnd and forkStart < forkEnd irrespective of alignOrn and forkDirn.
#             - pauseSite is the genomic coordinate of the pause site.
#             - paramStrNoPause is a string of the format high_low_offset_width, these are the parameters of the
#               no-cut sigmoid e.g.: 0.64_0.1_30934.820415496826_-3000. width < 0 for a left-moving fork and > 0 for
#               a right-moving fork.
#             - leftPieceXmax and rightPieceXmin correspond to coordinates of the ends of the cut pieces of the fork
#               aligned to the sigmoid measured relative to the inflection point of the sigmoid. In a right-moving fork,
#               leftPieceXmax is the distance of the pause site from the inflection point and for a left-moving fork,
#               -1 * rightPieceXmin is this quantity (with a negative sign as the reference curve is now reflected
#               about the y-axis).
#             - width is the width of the sigmoid (< 0 for a left-moving fork and > 0 for a right-moving fork).
#             - any column whose name starts with "keep" must have the values True or False.
#               These columns are used to filter out rows from the pause file based on various criteria.
# use_only_these_keep_columns: comma-separated list of columns to use for filtering the pause file for distributing
#                              probabilities. e.g.: "keep_fs,keep_fl" (with or without quotes).
#                              This is a subset of the columns whose name starts with "keep" in the pause file.
#                              This criterion is only used while assigning probabilities, not while calculating them.
#                              The reason we need this is because while assigning probabilities, not all keep* columns
#                              are relevant.
#                              - For example, we may have marked a fork for rejection because the pauseSite
#                                was too close to one end. But, probabilities can still be calculated for sites which
#                                are not too close to the end. So, this specific keep column has to be ignored while
#                                assigning probabilities and hence, should not be used for filtering the pause file and
#                                should not be included in this list.
#                              - Example 2: a keep column may mark a fork for discard because it is too short.
#                                This criterion should stop us from calculating probabilities for the fork.
#                                So, this specific keep column should be used for filtering the pause file and should
#                                be included in this list.
# fasta: path to fasta file containing the reference genome. Candidate pause sites are thymidines and we can only
#        know where those are using the fasta file. NOTE: thymidines are hard-coded somewhere in this script or in
#        an associated script, so if you want to use a different base, you'll have to change that.
# decimate_T_optional: only output probability vs fork coordinate every decimate_T_optional thymidines.
#                      default is 20. saves compute time.

# logic
# -----
# - find all pauses (i.e. all rows with all keep* = True irrespective of the value of use_only_these_keep_columns)
# - calculate the distance of each pause from the inflection point and construct a p.d.f. of these distances
# - keep all forks without a pause using use_only_these_keep_columns
# - assign a probability to each site along those forks using the p.d.f. of distances from the inflection point
# - reject pause sites suitably using our filtration criteria in process_model_pause_results_analysis.py
#   (NOTE: not all criteria in that script can be applied here e.g. we don't have a pause duration here, so applying a
#    duration-based criterion does not make sense. We can turn on/off different parts of that script by ensuring
#    that some columns are not available in our probability file).
# - convert the un rejected pause sites to bed format as explained in the output section below.

# outputs
# -------
# Bed file sent to standard output with the format:
# - Eight tab-separated columns with comments starting with '#' and no column names.
# - First six columns are the standard bed columns of contig, start, end, name, score, strand and names are read ids.
# - Score is always set to 1000. If you see a different value, it is a bug.
# - The seventh column is numerical and is the pause sensitivity at that base at that read.
# - The eighth column is L or R depending on whether the current fork is left or right.
# - As there can be multiple forks per read, all information pertaining to a read may not appear contiguously.
# - Although the bed file entries are all single base, data is not output at every base along the read because
#   pauses can only be detected at thymidines and according to the decimate_T parameter.
#   Sigmoid profiles vary over ~kb scales, so dropping data on the ~10bp scale should not affect the results.
#   So we can save compute time.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -python -miller -jq

# load configuration
source config.sh

# set temporary directory and make it
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
keep_columns=${2:-}
fasta=${3:-}
decimate_T=${4:-20}

# check that the correct number of arguments were provided
if [ "$#" -lt 3 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: sbatch $0 pause_file keep_columns fasta decimate_T_optional"
    >&2 echo "pause_file: path to pause file obtained by the cut-and-align procedure using a reference sigmoid."
    >&2 echo "keep_columns: comma-separated list of columns to use for filtering the pause file."
    >&2 echo "fasta: path to fasta file containing the reference genome."
    >&2 echo "decimate_T_optional: only output probability vs fork coordinate every decimate_T_optional thymidines. "
    >&2 echo "                     default is 20. saves compute time."
    >&2 echo "For more details on what these parameters mean, see the comments in the script."
    >&2 echo "optional means that the argument is optional."
    exit 1;
fi

# check that the pause file exists
if [ ! -f "$pause_file" ]; then
    >&2 echo "ERROR: pause file $pause_file does not exist"
    exit 1;
fi

# check that the fasta file exists
if [ ! -f "$fasta" ]; then
    >&2 echo "ERROR: fasta file $fasta does not exist"
    exit 1;
fi

# check that decimate_T is an integer
if ! [[ "$decimate_T" =~ ^[0-9]+$ ]]; then
    >&2 echo "ERROR: decimate_T must be an integer"
    exit 1;
fi

# check that decimate_T is above 0
if [ "$decimate_T" -le 0 ]; then
    >&2 echo "ERROR: decimate_T must be above 0"
    exit 1;
fi

# make a list of required columns
# NOTE: detectIndex, pauseSite are automatically included in the list of required columns in the checking script
req_columns="$keep_columns,leftPieceXmax,rightPieceXmin,width"

# check that the input pause file is valid
if [ ! "$(< "$pause_file" python validate_pause_format.py --required-columns "$req_columns" )" == "valid" ]; then
  >&2 echo "Error: pause file is not valid."
  exit 1;
fi

# Step 1: get mean and s.d. of distances of pause sites from inflection point
# ===========================================================================
# logic is: first filter out all rows with any keep* = False, then calculate stats

tmp_stats_file=$(mktemp "$tmpDir"/pause_file_inflection_point.XXXXXXXXXX)
tmp_rel_pause_loc_list_file=$(mktemp "$tmpDir"/rel_pause_loc_list.XXXXXXXXXX)

# create a list of relative pause locations and send to file
# shellcheck disable=SC2016
# shellcheck disable=SC1010
mlr --tsv --skip-comments\
    put -e '$any_keep_false = 0;
            for (k, v in $*){
                if(k =~ "^keep"){
                    if(v == "False"){
                        $any_keep_false = 1;
                    }
                }
            }'\
    then filter '$any_keep_false == 0' then\
    put 'if($width > 0){
           $rel_pause_loc = $leftPieceXmax;
         } else {
           $rel_pause_loc = -$rightPieceXmin;
         }' then cut -f rel_pause_loc "$pause_file" > "$tmp_rel_pause_loc_list_file"

# get mean and sd of relative pause locations
mlr --itsv --ojson --skip-comments stats1 -a mean,stddev -f rel_pause_loc "$tmp_rel_pause_loc_list_file" >\
  "$tmp_stats_file"

pause_inflection_mean=$(jq -r '.[].rel_pause_loc_mean' "$tmp_stats_file")
pause_inflection_sd=$(jq -r '.[].rel_pause_loc_stddev' "$tmp_stats_file")

# do a sanity check that the sd is at least 1000 bp
if [ "$(echo "$pause_inflection_sd < 1000" | bc)" -eq 1 ]; then
  >&2 echo "ERROR: pause inflection point standard deviation is less than 1000"
  >&2 echo "This is probably an error. Please check the pause file and the script."
  exit 1;
fi

# Step 2: redistribute pause sites among non-paused forks using the p.d.f. of distances from the inflection point
# ===============================================================================================================
# logic is: first filter out all rows where any of the user-requested columns are False, then assign probabilities,
#           then filter using our model analysis script, and then convert to bed those positions where all keep* = True

# count total number of pauses to use in a normalization step later
pause_count=$(< "$pause_file" python convert_pause_file_to_bed.py --discardUnkeptPauses |\
                  grep -E -v '^browser|^track|^#' |\
                  mlr --tsv --implicit-csv-header --headerless-csv-output count)

# make a temp file to store non normalized bed data
tmp_non_normalized_bed_file=$(mktemp "$tmpDir"/non_normalized_bed.XXXXXXXXXX)

# make a temp file to store the output of the first step
tmp_prob_file=$(mktemp "$tmpDir"/prob_file.XXXXXXXXXX)

# convert keep_columns to a regex
keep_columns_regex=$(echo "($keep_columns)" | tr ',' '|')

# shellcheck disable=SC2016
# shellcheck disable=SC1010
< "$pause_file"\
  mlr --tsv --skip-comments put -s keep_columns_regex="$keep_columns_regex" -e\
           '$given_keep_false = 0;
            for (k, v in $*){
                if(k =~ "@keep_columns_regex"){
                    if(v == "False"){
                        $given_keep_false = 1;
                    }
                }
            }' then filter '$given_keep_false == 0' |\
    python calculate_probabilities_using_sgm_pause_to_inflection_pdf_along_non_paused_forks.py\
      --gaussian-mean "$pause_inflection_mean" --gaussian-sd "$pause_inflection_sd" --fasta "$fasta"\
      --pause-to-inflection-list "$tmp_rel_pause_loc_list_file" --n-bins 20\
      --check-mean-sd-pdf-against-gaussian-with-tol 100 |\
        tee "$tmp_prob_file" |\
          mlr --tsv --skip-comments decimate -n "$decimate_T" |\
          python process_model_pause_results_analysis.py --keep_cols_requested keep_ep,keep_lp |\
            python convert_pause_file_to_bed.py --discardUnkeptPauses |\
              grep -v "^#" |\
                awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $9, $7 ~ /_L_/ ? "L" : "R"}'\
                   >"$tmp_non_normalized_bed_file"

# sum the seventh column
sum_seventh_column=$(awk 'BEGIN{sum = 0}{sum += $7}END{print sum}' "$tmp_non_normalized_bed_file")

# print comments output by the first step
head -n 100 "$tmp_prob_file" | grep "^#"

# normalize the seventh column and print output
awk -v sum_seventh_column="$sum_seventh_column" -v pause_count="$pause_count"\
  'BEGIN{OFS="\t"}{$7 = $7 * pause_count / sum_seventh_column; print}' "$tmp_non_normalized_bed_file" |\
    insert_calling_script_header "$pause_file" "$keep_columns" "$fasta"

# last step: clean up
# ===================
rm -rf "$tmpDir"