#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J regPauseStats
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# goal
# -----
# Measure various pause statistics for a region of interest and output them as a json object

# usage
#------
# bash calculate_region_pause_statistics.sh [-o output_bed_stats] op_dir prefix mod_bam pause_file bed_file relative_direction_option\
#     genome_size_bp_optional pause_sensitivity_bedgraphs_prefix_optional restrict_fork_direction_optional\
#     delete_new_mod_bam_made_on_the_fly_optional AT_ratio_optional fasta_optional
# NOTE: can use sbatch in place of bash
# NOTE: [] denotes optional arguments. If you don't want to specify it, ignore stuff between [].
#       If you want to specify it, remove the [] and write the text within the [].
# -o output_bed_stats: (default /dev/null) write some statistics about the input bed file (containing the ROI) to this
#                      file.
# op_dir: output directory
# prefix: prefix for output files
# mod_bam: modified bam file containing all reads and their modification information
# pause_file: pause file containing all forks and their pause information
# bed_file: bed file containing regions of interest
# relative_direction_option: (default "all") type of relative directions between fork and region to allow,
#                            set to "all", "co-directional" or "head-on".
# genome_size_bp_optional: (default 0) size of the genome in bps, used to calculate fraction of genome covered by
#                          bed file.
#                          NOTE: by default, this calculation is not performed.
#                          NOTE: we don't automatically do this calculation using a fasta file as there are assumptions
#                                like how many rDNA regions to use, exclude chrM etc.
# pause_sensitivity_bedgraphs_prefix_optional: (default "" i.e. unused) prefix for pause sensitivity bedgraphs.
#                                               Let's say this prefix is set to /path/to/file, then the script
#                                               will look for files like: /path/to/file.all.bedgraph,
#                                               /path/to/file.left.bedgraph, /path/to/file.right.bedgraph.
#                                               These files are four-column text-separated files with no column names,
#                                               comments starting with '#' and the columns: contig, start, end, value.
#                                               The value column is the pause sensitivity value for that position.
#                                               See run_get_sgm_pause_sensitivities.sh for how to generate such files
#                                               or other details. By default, this calculation is not performed.
# restrict_fork_direction_optional: (default "" i.e. no restriction) options are: "L", "R","lead","lag","" or unset.
#                                    if set to "L" or "R", restrict all calculations to only left- or right-moving
#                                    forks respectively. If set to "lead" or "lag", restrict all calculations to only
#                                    forks corresponding to the leading or lagging strand synthesis respectively.
#                                    (In reality, forks perform both leading and lagging strand synthesis; each
#                                     fork synthesizes sections belonging to two different molecules.
#                                     As our data is all single molecule, it makes sense to classify forks as performing
#                                     leading or lagging strand synthesis based on the direction of the fork movement
#                                     and the alignment of the read to the reference genome.).
# delete_new_mod_bam_made_on_the_fly_optional: (default "FALSE") if set to "TRUE", delete the modified bam file made on
#                                              the fly which was made by intersecting the mod bam and bed file.
#                                              This speeds up the job and saves disk space.
# AT_ratio_optional: (default 0) ratio of number of A and T bases to total number of bases in the reference genome.
#                    By default, this option is unused.
#                    How this is used: when we calculate null hypotheses like "if pauses occur uniformly across
#                    the genome, how many pauses would we expect in the region of interest?", we need to account for
#                    our method calling pauses only at thymidines. So in regions of high AT content, we would expect
#                    more pauses (I say AT instead of just T because pauses can happen on the reference strand or
#                    its complement and we usually include pauses on both strands in our intersection calculations).
#                    So we report something like "Fraction of the AT genome covered by the ROI" in the output json.
#                    This quantity is not reported if the AT_ratio is not set.
#                    We do not calculate this from the fasta file as there are repetitive regions like rDNA
#                    that are not in the file, and we may not want to include regions like the mitochondrial genome
#                    in the calculation etc. etc.
#                    So, we expect the user to input this quantity after having calculated it themselves.
# fasta_optional: (default "") fasta file to calculate base content to calculate ratio of AT bases to total bases in the
#                 region of interest in the reference genome. By default, this calculation is not performed.

# outputs
# -------
# produces a json object, an array which looks like this:
#[
#{
#    "value": 1000,
#    "desc": "blah_blah_blah",
#    "script_and_args": "calculate_region_pause_statistics.sh blah blah",
#    "commit": "fa8632e",
#    "date": "Mon 27 Mar 17:20:02 BST 2023",
#    "author": "Blah blah <blah@blah.com>",
#    "latex_desc": "Blah blah blah"
#},...
#]

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages
source load_package.sh -miller -python -samtools -jq > /dev/null 2>&1

# load configuration
source config.sh

# make a temp directory in scratch
mkdir -p "${config[scratchDir]:-}"/tmp

# process flags
output_bed_stats=/dev/null
while getopts ":o:" opt; do
  case $opt in
    o) output_bed_stats="$OPTARG"
    ;;
    \?) : # do nothing
    ;;
  esac
done

# shift arguments to exclude parsed options
shift $((OPTIND-1))

# check if there are enough arguments
if [ "$#" -lt 5 ]; then
  >&2 echo "Error: incorrect number of arguments"
  >&2 echo "Usage: bash $0 [-o output_bed_stats] op_dir prefix mod_bam pause_file bed_file relative_direction_optional \ "
  >&2 echo "       genome_size_bp_optional pause_sensitivity_bedgraphs_prefix_optional \ "
  >&2 echo "       restrict_fork_direction_optional delete_new_mod_bam_made_on_the_fly_optional AT_ratio_optional \ "
  >&2 echo "       fasta_optional"
  >&2 echo "NOTE: \ means the command is continued on the next line for readability"
  >&2 echo "NOTE: [] denotes optional arguments. If you don't want to specify it, ignore stuff between []."
  >&2 echo "      If you want to specify it, remove the [] and write the text within the []."
  >&2 echo " -o output_bed_stats: (default /dev/null) write some statistics about the input bed file (containing the ROI) to the file at output_bed_stats"
  >&2 echo "op_dir: output directory"
  >&2 echo "prefix: prefix for output files"
  >&2 echo "mod_bam: modified bam file containing all reads and their modification information"
  >&2 echo "pause_file: pause file containing all forks and their pause information"
  >&2 echo "bed_file: bed file containing regions of interest"
  >&2 echo "relative_direction_option: (default \"all\") type of relative directions between fork and region to allow, set to \"all\", \"co-directional\" or \"head-on\"."
  >&2 echo "genome_size_bp_optional: (default 0 i.e. not used) size of the genome in bps "
  >&2 echo "                         (remember to include many rDNAs and exclude mitochondrial genome as need be),"
  >&2 echo "                         used to calculate fraction of genome covered by bed file."
  >&2 echo "pause_sensitivity_bedgraphs_prefix_optional: (default \"\" i.e. not used) prefix for pause sensitivity bedgraphs."
  >&2 echo "                                            Let's say this prefix is set to /path/to/file, then the script"
  >&2 echo "                                            will look for files like: /path/to/file.all.bedgraph,"
  >&2 echo "                                            /path/to/file.left.bedgraph, /path/to/file.right.bedgraph."
  >&2 echo "                                            For more details, see the script."
  >&2 echo "restrict_fork_direction_optional: (default \"\" i.e. no restriction) set to \"L\" or \"R\" or \"lead\" or \"lag\","
  >&2 echo "                                    if set to \"L\" or \"R\", restrict all calculations to only left- or right-moving"
  >&2 echo "                                    forks respectively. If set to \"lead\" or \"lag\", restrict all calculations to only"
  >&2 echo "                                    forks corresponding to the leading or lagging strand synthesis respectively."
  >&2 echo "                                    As our data is single-molecule, forks _can_ be classified as leading or lagging"
  >&2 echo "delete_new_mod_bam_made_on_the_fly_optional: (default \"FALSE\") if set to \"TRUE\", delete the "
  >&2 echo "                                             modified bam file made on the fly which was made by "
  >&2 echo "                                             intersecting the mod bam and bed file. "
  >&2 echo "                                             This speeds up the job and saves disk space."
  >&2 echo "AT_ratio_optional: (default 0 i.e. unused) ratio of number of A and T bases to total number of bases in the"
  >&2 echo "                   reference genome."
  >&2 echo "fasta_optional: (default \"\") fasta file to calculate base content to see if the region is AT-rich etc."
  >&2 echo "               By default, this calculation is not performed."
  exit 1
fi

# load variables from the command line
op_dir=$1
prefix=$2
mod_bam=$3
pause_file=$4
bed_file=$5
relative_direction_option=${6:-"all"}
genome_size_bp=${7:-0}
pause_sens_file_prefix=${8:-""}
restrict_fork_direction=${9:-""}
delete_new_mod_bam_made_on_the_fly=${10:-"FALSE"}
AT_ratio=${11:-0}
fasta_file=${12:-/dev/null}

# check that mod_bam, pause_file, and bed_file exist
if [ ! -f "$mod_bam" ] || [ ! -f "$pause_file" ] || [ ! -f "$bed_file" ]; then
  >&2 echo "Error: one or more input files do not exist"
  exit 1
fi

# check that mod_bam is indexed
if [ ! -f "$mod_bam".bai ]; then
  >&2 echo "Error: mod_bam is not indexed"
  exit 1
fi

# if pause_sens_file_prefix is not empty, check that the three bedgraph files exist
# if they do exist, set flag to perform pause sensitivity calculations
is_pause_sens_calc=false
if [ -n "$pause_sens_file_prefix" ]; then

  # we check for three files here: all, left, right
  # if the user wants to restrict the fork direction, then we check for four additional files:
  #   left.minus, left.plus, right.plus, right.minus
  if [ "$restrict_fork_direction" == "lead" ] || [ "$restrict_fork_direction" == "lag" ]; then
    suffix_list=(left.minus left.plus right.plus right.minus all left right)
  else
    suffix_list=(all left right)
  fi

  for suffix in "${suffix_list[@]}"; do

    sensitivity_file="$pause_sens_file_prefix"."$suffix".bedgraph

    if [ ! -f "$sensitivity_file" ]; then
      >&2 echo "Error: pause sensitivity bedgraph file $pause_sens_file_prefix.$suffix.bedgraph does not exist";
      exit 1;
    elif [ -f "$fasta_file".fai ]; then
      validity=$(< "$sensitivity_file" grep -v -E '^browser|^track|^#' | awk 'BEGIN{OFS="\t"}{print $1, $2, $3}' |\
        python validate_bed_against_fai.py --check-genome-cov-to-tol 100  --check-no-overlap "$fasta_file".fai)

      if [ "$validity" != "valid" ]; then
        >&2 echo "ERROR: sensitivity file $sensitivity_file is not valid"
        >&2 echo "This may not mean that the file is an invalid bedgraph, but it does not meet the requirements "
        >&2 echo "for this script."
        exit 1;
      fi
    fi
  done

  # set flag to perform pause sensitivity calculations if the checks above pass
  is_pause_sens_calc=true

fi

# check that the input bed file is valid
if [ ! "$(< "$bed_file" python validate_bed_format.py --no-dot-strand --exactly-six-columns --allow-float-score)" == "valid" ]; then
  >&2 echo "Error: bed file is not valid."
  exit 1;
fi

# if the fasta file exists, then validate the bed file against the fasta .fai file
if [ -f "$fasta_file" ]; then

  if [ ! -f "$fasta_file".fai ]; then
    >&2 echo "Error: fasta index file $fasta_file.fai does not exist."
    exit 1;
  fi

  # bed file must have contigs only in the fasta index file
  if [ ! "$(< "$bed_file" python validate_bed_against_fai.py "$fasta_file".fai)" == "valid" ]; then
    >&2 echo "Error: bed file must have contigs only in the fasta index file."
    exit 1;
  fi

fi

# check that the input pause file is valid
if [ ! "$(< "$pause_file" python validate_pause_format.py)" == "valid" ]; then
  >&2 echo "Error: pause file is not valid."
  exit 1;
fi

# create output directory if it does not exist
mkdir -p "$op_dir"
op_dir=$(cd "$op_dir" || exit; pwd)

# if a fork-direction restriction is requested, make a new pause file with only forks moving in the requested direction,
# and set the pause file to this new file
if [ "$restrict_fork_direction" == "L" ] || [ "$restrict_fork_direction" == "R" ]; then

  tmp_file_forks_restricted=$(mktemp -p "${config[scratchDir]:-}"/tmp)
  # shellcheck disable=SC2016
  mlr --tsv --skip-comments filter '$detectIndex =~ "_'"$restrict_fork_direction"'_"' "$pause_file" > \
    "$tmp_file_forks_restricted"
  pause_file="$tmp_file_forks_restricted"

elif [ "$restrict_fork_direction" == "lead" ]; then

  tmp_file_forks_restricted=$(mktemp -p "${config[scratchDir]:-}"/tmp)
  # FLAG: LEAD LAG DISTINCTION
  # shellcheck disable=SC2016
  mlr --tsv --skip-comments filter '$detectIndex =~ "_fwd_R_" || $detectIndex =~ "_rev_L_"' "$pause_file" > \
    "$tmp_file_forks_restricted"
  pause_file="$tmp_file_forks_restricted"

elif [ "$restrict_fork_direction" == "lag" ]; then

  tmp_file_forks_restricted=$(mktemp -p "${config[scratchDir]:-}"/tmp)
  # FLAG: LEAD LAG DISTINCTION
  # shellcheck disable=SC2016
  mlr --tsv --skip-comments filter '$detectIndex =~ "_fwd_L_" || $detectIndex =~ "_rev_R_"' "$pause_file" > \
    "$tmp_file_forks_restricted"
  pause_file="$tmp_file_forks_restricted"

else

  restrict_fork_direction=""

fi

SCRIPT_AND_ARGS="$0 $*"

# set output mod bam file
op_mod_bam="$op_dir"/"$prefix".mod.bam

# set output pause file
op_all_pause_file="$op_dir"/"$prefix"_all_pauses # here, "all" means pauses whether we have confidence in them or not
op_no_pause_file="$op_dir"/"$prefix"_no_pauses   # forks with no pauses in them
op_elsewhere_pause_file="$op_dir"/"$prefix"_pauses_elsewhere # forks that pass through region but pause elsewhere
op_region_pause_file="$op_dir"/"$prefix"_pauses_in_region    # forks that pass through region and pause in it

# function that outputs json objects
function create_json() {
    # if $1 is a number, then it is a value, otherwise it is a string
    if [[ $1 =~ ^[0-9.-]+$ ]]; then
        value=$1
    else
        value="\"$1\""
    fi
    desc="$2"
    latex_desc="$3"
    script_and_args="$SCRIPT_AND_ARGS"
    commit=${COMMITSTR:-NA}
    date=${TIMENOW:-NA}
    author="${config[name]:-NA} <${config[email]:-NA}>"

    cat <<EOF
{
    "value": $value,
    "desc": "$desc",
    "script_and_args": "$script_and_args",
    "commit": "$commit",
    "date": "$date",
    "author": "$author",
    "latex_desc": "$latex_desc"
},
EOF
}

# function that counts pause entries
report_count_pause_file() {
    local latex_desc="$3"
    local desc="$2"
    create_json "$(mlr --itsv --otsv --skip-comments  --headerless-csv-output count "$1")" "$desc" "$latex_desc"
}


# start output of json object
echo "["

create_json "$(samtools view -c "$mod_bam")" "${prefix}_n_reads" "Total number of reads in mod bam file"

report_count_pause_file "$pause_file" "${prefix}_n_forks" "Total number of forks in pause file"

# create a temporary file with only the true pauses by discarding any pause where any keep column is not True
tmp_true_pauses_file=$(mktemp -p "${config[scratchDir]:-}"/tmp)
{
  grep -m 1 -v '^#' "$pause_file";

< "$pause_file" python convert_pause_file_to_bed.py --discardUnkeptPauses |\
    grep -E -v '^browser|^track|^#' | cut -f7-

} > "$tmp_true_pauses_file"

# shellcheck disable=SC2016
# shellcheck disable=SC1010
n_total_pauses=$(< "$tmp_true_pauses_file" \
                  mlr --itsv --otsv --skip-comments --headerless-csv-output count)
create_json "$n_total_pauses" "${prefix}_n_pauses" "Total number of pauses in pause file $ N_p^{tot} = $"

# first make the mod bam file containing only reads that overlap the ROI and count the reads
if [ "$delete_new_mod_bam_made_on_the_fly" == "TRUE" ]; then
  read_count="$(samtools view -c "$mod_bam" --regions-file "$bed_file")"
else
  samtools view -b "$mod_bam" --regions-file "$bed_file" > "$op_mod_bam".temp
  samtools sort -o "$op_mod_bam" -T "$op_dir" "$op_mod_bam".temp
  samtools index "$op_mod_bam"
  rm "$op_mod_bam".temp
  read_count="$(samtools view -c "$op_mod_bam")"
fi

create_json "$read_count" "${prefix}_n_reads_roi"\
    "Total number of reads in mod bam file that overlap the ROI"

# intersect the pause file with the regions of interest
if [ "$relative_direction_option" == "all" ]; then
  bash intersect_pauses_with_bed.sh "$bed_file" "$pause_file" "$op_all_pause_file" "$op_no_pause_file"\
    "$op_elsewhere_pause_file" "$op_region_pause_file"
elif [ "$relative_direction_option" == "co-directional" ]; then
  bash intersect_pauses_with_bed.sh -fs "$bed_file" "$pause_file" "$op_all_pause_file" "$op_no_pause_file"\
    "$op_elsewhere_pause_file" "$op_region_pause_file"
elif [ "$relative_direction_option" == "head-on" ]; then
  bash intersect_pauses_with_bed.sh -fS "$bed_file" "$pause_file" "$op_all_pause_file" "$op_no_pause_file"\
    "$op_elsewhere_pause_file" "$op_region_pause_file"
else
  >&2 echo "Error: relative_direction_option must be set to 'all', 'co-directional' or 'head-on'"
  exit 1;
fi

report_count_pause_file "$op_all_pause_file" "${prefix}_n_forks_roi"\
  "Total number of forks in pause file that overlap the ROI"

report_count_pause_file "$op_no_pause_file" "${prefix}_n_forks_roi_no_pause"\
  "Total number of forks in pause file that overlap the ROI but have no pauses"

report_count_pause_file "$op_elsewhere_pause_file" "${prefix}_n_pauses_roi_elsewhere"\
  "Total number of pauses in pause file occurring on forks that overlap the ROI but not in the ROI"

latex_desc="Total number of pauses in pause file occurring on forks that overlap the ROI and in the ROI $ N_p^{reg} = $"
latex_desc="${latex_desc}NOTE: In this counting, pauses are not counted multiple times. For example, let's say a pause"
latex_desc="${latex_desc}      occurs in a region that is specified twice in the ROI bed file."
latex_desc="${latex_desc}      In spite of this, the pause is only counted once in this counting."

n_pause_region=$(mlr --tsv --skip-comments  --headerless-csv-output count "$op_region_pause_file")
create_json "$n_pause_region" "${prefix}_n_pauses_roi" "$latex_desc"

# If requested, perform pause sensitivity aggregation calculations.
# There are 15 cases here: 3 relative direction options x 5 fork direction options, and we use
# standard tools in the bedtools suite to do our calculations in the calculate_bed_overlap_bedgraph_signal_stats.sh
# script.
# The 3 relative direction options are: all, co-directional, head-on.
# The 5 fork direction options are: L, R, lead (which means R+ & L-), lag (which means R- & L+), "" (i.e. none).
# To leverage tools from the bedtools suite, we assign any L bedgraphs (just L, L+ or L-) to the minus strand
# and similarly for R bedgraphs. In this way, any co-directional intersect is a same-strand intersect and any
# head-on intersect is an opposite-strand intersect.

# FLAG: LEAD LAG DISTINCTION
if [ "$is_pause_sens_calc" == "true" ]; then

  # Set the relative strand (same/opposite) and the strand (p/m) for the bedgraph calculations
  # We act like left fork data is on the minus strand, and right fork data is on the plus strand
  # so that we can use same-strand or opposite-strand intersects for co-directional and head-on respectively.
  strand_left=""
  strand_right=""
  if [ "$relative_direction_option" == "all" ]; then
    : # do nothing
  elif [ "$relative_direction_option" == "co-directional" ]; then
    strand_left="-s -m"
    strand_right="-s -p"
  elif [ "$relative_direction_option" == "head-on" ]; then
    strand_left="-S -m"
    strand_right="-S -p"
  else
    >&2 echo "Error: relative_direction_option must be set to 'all', 'co-directional' or 'head-on'"
    exit 1;
  fi

  # Set suffix for bedgraphs based on the fork direction restriction
  left_suffix="left"
  right_suffix="right"
  if [ "$restrict_fork_direction" == "lead" ]; then
    left_suffix="left.minus"
    right_suffix="right.plus"
  elif [ "$restrict_fork_direction" == "lag" ]; then
    left_suffix="left.plus"
    right_suffix="right.minus"
  elif [ "$restrict_fork_direction" == "L" ] || [ "$restrict_fork_direction" == "R" ] ||\
    [ -z "$restrict_fork_direction" ]; then
    : # do nothing
  else
    >&2 echo "Error: restrict_fork_direction must be set to 'L', 'R', 'lead', 'lag' or ''"
    exit 1;
  fi

  # perform calculations
  total_sens_overlap_L=0
  total_sens_overlap_R=0
  if [ ! "$restrict_fork_direction" == "R" ]; then
    # shellcheck disable=SC2086
    total_sens_overlap_L=$(bash calculate_bed_overlap_bedgraph_signal_stats.sh $strand_left "$bed_file"  \
      "$pause_sens_file_prefix"."$left_suffix".bedgraph | jq -r '.overlap_sum')
  fi
  if [ ! "$restrict_fork_direction" == "L" ]; then
    # shellcheck disable=SC2086
    total_sens_overlap_R=$(bash calculate_bed_overlap_bedgraph_signal_stats.sh $strand_right "$bed_file" \
      "$pause_sens_file_prefix"."$right_suffix".bedgraph | jq -r '.overlap_sum')
  fi
  total_sens_overlap=$(echo "scale=8; $total_sens_overlap_L + $total_sens_overlap_R" | bc)

  # report sensitivity in a json object
  latex_desc="Total number of pauses in ROI expected from sensitivity analysis $ N_p^{sens} = $"
  create_json "$total_sens_overlap" "${prefix}_n_pauses_roi_sens" "$latex_desc"

fi

create_json "$(< "$bed_file" grep -E -v '^browser|^track|^#' | awk 'BEGIN{a=0}{a+=$3-$2}END{print a}')"\
  "${prefix}_total_roi_length"\
  "Total length of ROI"

create_json "$(< "$bed_file" grep -c -E -v '^browser|^track|^#')"\
  "${prefix}_total_n_roi"\
  "Total number of regions in ROI"

# get some statistics associated with the bed file
tmp_bed_stats=$(mktemp -p "${config[scratchDir]:-}"/tmp)
bash bed_region_statistics.sh "$bed_file" "" "$fasta_file" | tee "$output_bed_stats" > "$tmp_bed_stats"

# get fork coordinates in bed format
tmp_file_forks=$(mktemp -p "${config[scratchDir]:-}"/tmp)
< "$pause_file" python convert_pause_file_to_bed.py --outputForks --LRtoPlusMinus > "$tmp_file_forks"
# perform calculation
if [ "$relative_direction_option" == "all" ]; then
  total_fork_overlap_length=$(bash calculate_bed_total_masked_length.sh "$bed_file" "$tmp_file_forks")
  total_merged_roi_length=$(jq -r '.total_length_bp_merge_ignore_strand' "$tmp_bed_stats")
elif [ "$relative_direction_option" == "co-directional" ]; then
  total_fork_overlap_length=$(bash calculate_bed_total_masked_length.sh -s "$bed_file" "$tmp_file_forks")
  total_merged_roi_length=$(jq -r '.total_length_bp_merge_same_strand' "$tmp_bed_stats")
elif [ "$relative_direction_option" == "head-on" ]; then
  total_fork_overlap_length=$(bash calculate_bed_total_masked_length.sh -S "$bed_file" "$tmp_file_forks")
  total_merged_roi_length=$(jq -r '.total_length_bp_merge_same_strand' "$tmp_bed_stats")
else
  >&2 echo "Error: relative_direction_option must be set to 'all', 'co-directional' or 'head-on'"
  exit 1;
fi

latex_desc="Total length of ROI after merging overlapping regions. $ L_{ROI,merge} = $"
latex_desc="${latex_desc}NOTE: this could be slightly different depending on whether "
latex_desc="${latex_desc}the relative direction option was set to all or otherwise."

create_json "$total_merged_roi_length"\
  "${prefix}_total_merged_roi_length"\
  "$latex_desc"

if [ "$genome_size_bp" -gt 0 ]; then

  create_json "$genome_size_bp"\
    "${prefix}_total_genome_length"\
    "Total length of genome $ L_{genome} = $"

  fraction_genome_covered_by_merged_roi="$(echo "scale=8; $total_merged_roi_length/$genome_size_bp" | bc)"
  create_json "$fraction_genome_covered_by_merged_roi"\
    "${prefix}_total_merged_roi_length_fraction_genome"\
    "Fraction of genome covered by (merged) ROI $ f_{ROI} = L_{ROI,merge} / L_{genome} = $"

  if [ "$relative_direction_option" == "all" ]; then
    fraction_genome_covered_by_merged_roi_halved_if_need_be="$fraction_genome_covered_by_merged_roi"
  elif [ "$relative_direction_option" == "co-directional" ] || [ "$relative_direction_option" == "head-on" ]; then
    fraction_genome_covered_by_merged_roi_halved_if_need_be="$(echo "scale=8; $fraction_genome_covered_by_merged_roi/2" | bc)"
  else
    >&2 echo "Error: relative_direction_option must be set to 'all', 'co-directional' or 'head-on'"
    exit 1;
  fi

  latex_desc="Fraction of genome covered by (merged) ROI halved (assuming equal co-directional and head-on forks) $ f_{ROI,2} = $"
  latex_desc="${latex_desc}. If relative direction option is 'all', then same as above"
  create_json "$fraction_genome_covered_by_merged_roi_halved_if_need_be"\
    "${prefix}_total_merged_roi_length_fraction_genome_halved_if_need_be"\
    "$latex_desc"

  create_json "$(echo "sqrt($fraction_genome_covered_by_merged_roi_halved_if_need_be / $n_total_pauses)" | bc -l | sed 's/^\./0./')"\
    "${prefix}_total_merged_roi_length_fraction_genome_halved_if_need_be_sd"\
    "Expected error due to counting in the fraction above, $ \\\\sqrt{\\\\frac{f_{ROI,2}}{N_p^{tot}}} = $"

  # calculate the fraction of the AT genome covered by the ROI if the AT ratio is set and the fasta file exists
  if [ "$AT_ratio" != 0 ] && [ -f "$fasta_file" ]; then

    create_json "$AT_ratio"\
      "${prefix}_total_genome_fraction_AT"\
      "Fraction of AT bases in the genome (user-supplied) $ b_{AT,genome} = $ "

    T_num_same_strand="$(jq -r '.total_T_merge_same_strand' "$tmp_bed_stats")"
    T_num_opposite_strand="$(jq -r '.total_A_merge_same_strand' "$tmp_bed_stats")"
    T_num_both_strands_all_genome="$(echo "$AT_ratio * $genome_size_bp" | bc)"

    create_json "$T_num_same_strand"\
      "${prefix}_total_T_same_strand"\
      "Total number of T bases in (merged) ROI"

    create_json "$T_num_opposite_strand"\
      "${prefix}_total_T_opposite_strand"\
      "Total number of T bases in complementary strand of (merged) ROI"

    create_json "$T_num_both_strands_all_genome"\
      "${prefix}_total_T_both_strands_all_genome"\
      "Total number of T bases in both strands of genome"

    null_hypothesis_uniform_T_pause_all="$(echo "scale=8; ($T_num_same_strand + $T_num_opposite_strand)/$T_num_both_strands_all_genome" | bc)"
    null_hypothesis_uniform_T_pause_ho_or_cd_forks="$(echo "scale=8; $null_hypothesis_uniform_T_pause_all/2" | bc)"
    null_hypothesis_uniform_T_pause_ho_lead_or_cd_lag_forks="$(echo "scale=8; $T_num_opposite_strand/$T_num_both_strands_all_genome" | bc)"
    null_hypothesis_uniform_T_pause_ho_lag_or_cd_lead_forks="$(echo "scale=8; $T_num_same_strand/$T_num_both_strands_all_genome" | bc)"

    if [ "$relative_direction_option" == "all" ]; then
      null_T_hypothesis="$null_hypothesis_uniform_T_pause_all"
    elif { [ "$relative_direction_option" == "co-directional" ] && [ "$restrict_fork_direction" == "lead" ]; } || \
      { [ "$relative_direction_option" == "head-on" ] && [ "$restrict_fork_direction" == "lag" ]; }; then
      null_T_hypothesis="$null_hypothesis_uniform_T_pause_ho_lag_or_cd_lead_forks"
    elif { [ "$relative_direction_option" == "co-directional" ] && [ "$restrict_fork_direction" == "lag" ]; } || \
      { [ "$relative_direction_option" == "head-on" ] && [ "$restrict_fork_direction" == "lead" ]; }; then
      null_T_hypothesis="$null_hypothesis_uniform_T_pause_ho_lead_or_cd_lag_forks"
    elif [ "$relative_direction_option" == "head-on" ] || [ "$relative_direction_option" == "co-directional" ]; then
      null_T_hypothesis="$null_hypothesis_uniform_T_pause_ho_or_cd_forks"
    else
      >&2 echo "Error: relative_direction_option must be set to 'all', 'co-directional' or 'head-on'"
      exit 1;
    fi

    latex_desc="Expected pause fraction using null pause hypothesis of pauses at every available T $ f_{null,T} = $"
    latex_desc="${latex_desc}NOTE: number of available Ts depends on the relative direction option and the fork direction restriction."
    latex_desc="${latex_desc}For e.g.: if you choose head-on and leading strand synthesis, then a + element in the ROI will only see"
    latex_desc="${latex_desc}T bases on the - strand. For details, please see the calculations in the script."

    create_json "$null_T_hypothesis" \
      "${prefix}_null_hypothesis_for_uniform_T"\
      "$latex_desc"

    create_json "$(echo "sqrt($null_T_hypothesis / $n_total_pauses)" | bc -l | sed 's/^\./0./')"\
      "${prefix}_sd_null_hypothesis_for_uniform_T"\
      "s.d. of above null hypothesis $ \\\\sqrt{\\\\frac{f_{null,T}}{N_p^{tot}}} = $"

  fi

fi

latex_desc="Total length of sections of forks that overlap with the ROI $ N_f^{reg} = $. "
latex_desc="${latex_desc}NOTE: In this counting, forks are not counted multiple times. For example, let's say a fork"
latex_desc="${latex_desc}      overlaps with a region that is specified twice in the ROI bed file."
latex_desc="${latex_desc}      In spite of this, the length of the overlapping part of the fork is measured only once."

create_json "$total_fork_overlap_length"\
  "${prefix}_total_fork_overlap_length"\
  "$latex_desc"

total_fork_length=$(< "$tmp_file_forks" grep -v '^#' | awk 'BEGIN{a=0}{a+=$3-$2}END{print a}')
create_json "$total_fork_length"\
  "${prefix}_total_fork_length"\
  "Total length of forks in the pause file $ N_f^{tot} = $"

latex_desc="Ratio of total length of sections of forks that overlap the ROI to total length of forks $ N_f^{reg} / N_f^{tot} = $"
latex_desc="${latex_desc} NOTE: this is the expected ratio of pauses in ROI to total number of pauses if region is typical."
create_json "$(echo "$total_fork_overlap_length / $total_fork_length" | bc -l | sed 's/^\./0./')"\
  "${prefix}_ratio_fork_overlap_length_to_fork_length"\
  "$latex_desc"

create_json "$(echo "sqrt($total_fork_overlap_length / ($total_fork_length * $n_total_pauses))" | bc -l | sed 's/^\./0./')"\
  "${prefix}_ratio_fork_overlap_length_to_fork_length_th_sd"\
  "(theoretical) s.d. in expected ratio using counting statistics $ \\\\sqrt{\\\\frac{N_f^{reg}}{N_f^{tot} N_p^{tot}}} = $"

if [ "$is_pause_sens_calc" == "true" ]; then
  create_json "$(echo "scale=8; $total_sens_overlap / $n_total_pauses" | bc -l | sed 's/^\./0./')"\
    "${prefix}_ratio_sens_region_to_total_pauses"\
    "Ratio of number of pauses expected from sensitivity in ROI to total number of pauses $ N_p^{sens} / N_p^{tot} = $"

  create_json "$(echo "scale=8; sqrt($total_sens_overlap) / $n_total_pauses" | bc -l | sed 's/^\./0./')"\
    "${prefix}_ratio_sens_region_to_total_pauses_sd"\
    "s.d. in number of pauses expected from sensitivity in ROI to total number of pauses $ \\\\frac{\\\\sqrt{N_p^{sens}}}{N_p^{tot}} = $"
fi

create_json "$(echo "$n_pause_region / $n_total_pauses" | bc -l | sed 's/^\./0./')"\
  "${prefix}_ratio_pause_region_to_total_pauses"\
  "Ratio of number of pauses in ROI to total number of pauses $ N_p^{reg} / N_p^{tot} = $"

# Report statistics related to pause duration in region
tmp_duration_stats=$(mktemp -p "${config[scratchDir]:-}"/tmp)

mlr --itsv --ojson --skip-comments stats1 -a mean,stddev -f pauseDuration "$op_region_pause_file" \
  > "$tmp_duration_stats"

create_json "$(jq -r '.[].pauseDuration_mean // "NA"' "$tmp_duration_stats")"\
  "${prefix}_mean_duration"\
  "Mean duration of valid pauses in region (in kb)"

create_json "$(jq -r '.[].pauseDuration_stddev // "NA" | if . == "" then "NA" else . end' "$tmp_duration_stats")"\
  "${prefix}_sd_duration"\
  "S.d. of duration of valid pauses in region (in kb)"

# Report statistics related to pause duration in the whole pause file (restricted to valid pauses)
mlr --itsv --ojson --skip-comments stats1 -a mean,stddev -f pauseDuration "$tmp_true_pauses_file" \
  > "$tmp_duration_stats"

create_json "$(jq -r '.[].pauseDuration_mean // "NA"' "$tmp_duration_stats")"\
  "${prefix}_mean_duration_all_pauses"\
  "Mean duration of valid pauses in the entire pause file (kb)"

create_json "$(jq -r '.[].pauseDuration_stddev // "NA" | if . == "" then "NA" else . end' "$tmp_duration_stats")"\
  "${prefix}_sd_duration_all_pauses"\
  "S.d. of duration of valid pauses in the entire pause file (kb)"

# need an empty object at the end
echo "{}"

# close the json file
echo "]"

# remove temporary files
rm "$tmp_file_forks"
rm "$tmp_duration_stats"
rm "$tmp_bed_stats"
rm "$tmp_true_pauses_file"
if [ "$restrict_fork_direction" == "L" ] || [ "$restrict_fork_direction" == "R" ] ||\
  [ "$restrict_fork_direction" == "lead" ] || [ "$restrict_fork_direction" == "lag" ]; then
  rm "$tmp_file_forks_restricted"
fi