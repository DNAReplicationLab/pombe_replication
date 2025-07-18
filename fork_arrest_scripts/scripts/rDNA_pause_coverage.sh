#!/bin/bash

# goal
# -----
# Make a bunch of bedgraphs from pause files to see if there's a spatial pattern to pause occurrence

# usage
#------
# bash <script_name.sh> ref_fasta contig_name num_repeat_units pause_file bed_file_with_features output_dir nascent_bed\
#     base number pos_boundary parent_bed
# use bash instead of sbatch to run locally
# ref_fasta: path to reference fasta file of the rDNA repeat
# contig_name: name of contig of interest, usually something like rDNA_22_repeat
# num_repeat_units: number of repeat units in ref_fasta, in the example above it's 22
# pause_file: path to pause file (see details)
# bed_file_with_features: path to bed file with features (see details)
# output_dir: path to output directory, preferably this is a directory with nothing in it
# nascent_bed: (optional) bed file with alignment coordinates of reads identified as nascent (see details)
# base: (optional) base to create tilings on the genome, can be set to A, C, G, T or N (default N)
# number: (optional) number of bases to create tilings on the genome, can be set to any positive integer (default 13)
# pos_boundary: (optional) tile genome such that windows boundaries are forced at this position (default is
#                no such forcing is done and any suitable tiling of the genome is ok)
# parent_bed: (optional) bed file with alignment coordinates of reads identified as parent (see details)

# details
#---------
# pause_file is tab-delimited with headers with comments starting with '#' and column names.
# pause files must contain two columns: the detectIndex and the pauseSite.
#   detectIndex has the format readID_contig_start_end_orientation_direction_startFork_endFork.
#     orientation = fwd or rev, direction = L or R, startFork and endFork are the coordinates of the fork.
#     start always < end, startFork always < endFork both irrespective of orientation and direction.
#  pauseSite is the coordinate of the pause site, and must lie within startFork and endFork.

# bed_file_with_features is a bed file with features on the rDNA.
# an example is given below.
# NOTE: replace spaces with tabs in the example below.
# NOTE: do not assume the data below is real data, it is just an example.
# rDNA_22_repeat	-110	6747	RDN37-1	0	-
# rDNA_22_repeat	7991	8111	RDN5-1	0	+
# rDNA_22_repeat	7306	7412	ARS1200-1	0	-
# rDNA_22_repeat	8855	8874	RFB1	0	+
# rDNA_22_repeat	8915	8934	RFB2	0	+
# rDNA_22_repeat	8952	8972	RFB3	0	+

# nascent_bed is a bed file with the coordinates of reads identified as nascent reads.
# It is output by the file rDNA_detectSummary.py.
# The format is BED6 i.e. 6 tab-separated columns of contig, start, end, name, score, strand.
# The score could be a floating point number.

# parent_bed is same as above for reads identified as parental i.e. BrdU is too low to be considered nascent.

# outputs
# -------
# several bedgraphs and bed files are sent to output_dir, so its preferable to use an empty directory

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages
source load_package.sh -python

# load configuration
source config.sh

# assign arguments to variables
ref_fasta=$1
ref_fasta_index="$ref_fasta".fai
contig_name=$2
num_repeat_units=$3
pause_file=$4
bed_file_with_features=$5
output_dir=$6
nascent_bed=${7:-}
window_base=${8:-N}
window_number_bases=${9:-13}
window_force_boundary=${10:-NA}
parent_bed=${11:-}

# check that the correct number of arguments were provided
if [ "$#" -lt 6 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash $0 ref_fasta contig_name num_repeat_units pause_file bed_file_with_features output_dir \ "
    >&2 echo "               nascent_bed_optional base_optional number_optional pos_boundary_optional parent_bed_optional"
    >&2 echo "NOTE: The suffix _optional means that the parameter is optional."
    >&2 echo "NOTE: Optional parameters needn't be provided but if they are, they must be in the order shown above."
    >&2 echo "      What this usually means is, if you want to provide the 2nd optional argument, "
    >&2 echo "      you must provide the 1st optional argument too, but may skip the third, fourth,... etc."
    >&2 echo "NOTE: \  means the command is continued on the next line"
    >&2 echo "NOTE: A brief description of parameters is given below. For more details, see comments in the script."
    >&2 echo "ref_fasta is the path to the reference fasta file containing just the rDNA region repeated many times"
    >&2 echo "contig_name is the name of the contig of interest, usually something like rDNA_22_repeat"
    >&2 echo "num_repeat_units is the number of repeat units in ref_fasta, in the example above it's 22"
    >&2 echo "pause_file is the path to the pause file with the rDNA pauses in it"
    >&2 echo "bed_file_with_features is the path to the bed file with features on the rDNA"
    >&2 echo "output_dir is the path to the output directory, preferably this is a directory with nothing in it"
    >&2 echo "nascent_bed_optional is the path to the bed file with the coordinates of reads identified as nascent reads"
    >&2 echo "base_optional is the base to create tilings on the genome, can be set to A, C, G, T or N (default N)"
    >&2 echo "number_optional is the number of bases to create tilings on the genome, can be set to any positive integer (default 13)"
    >&2 echo "pos_boundary_optional is the position to tile the genome such that windows boundaries are forced at this position"
    >&2 echo "                      (default is no such forcing is done and any suitable tiling of the genome is ok)"
    >&2 echo "parent_bed_optional is the path to the bed file with the coordinates of reads identified as parental"
    exit 1;
fi

# make output directory if it does not exist
mkdir -p "$output_dir"
output_dir=$(realpath "$output_dir")

# check that the fasta and fasta index files exist
if [ ! -f "$ref_fasta_index" ] || [ ! -f "$ref_fasta" ]; then
    >&2 echo "ERROR: fasta index file and/or fasta file do not exist"
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

# check that the pause file exists
if [ ! -f "$pause_file" ]; then
    >&2 echo "ERROR: pause file does not exist"
    exit 1;
fi

# check that the pause file is valid
if [ ! "$(< "$pause_file" python validate_pause_format.py)" == "valid" ]; then
  >&2 echo "Error: pause file is not valid."
  exit 1;
fi

# if nascent and/or parental bed files are provided, check that they are valid
for file in "$nascent_bed" "$parent_bed"; do

  if [ -f "$file" ]; then

    if [ ! "$(< "$file" python validate_bed_format.py --six-columns --allow-float-score --require-uuid)" == "valid" ]; then
      >&2 echo "Error: $file is invalid."
      exit 1;
    fi

  fi

done

# set up the script calling information header in a function
# ==========================================================
insert_calling_script_header() {
  sed '1i'\
'# from commit '"${COMMITSTR:-NA}"' generated at '"${TIMENOW:-NA}"' by '"${config[name]:-NA}"' <'"${config[email]:-NA}"'>\n'\
'# script: '"$0"'\n'\
'# arguments: '"$*"'\n'\
"# slurm job name: ${SLURM_JOB_NAME:-NA}"
}

# get the length of the contig of interest and the length of each repeat unit
# ===========================================================================

contig_length=$(< "$ref_fasta_index" grep "$contig_name" | awk '{print $2}')
repeat_length=$((contig_length/num_repeat_units))
repeat_length_remainder=$((contig_length % num_repeat_units))

# check that the contig length is a multiple of the number of repeat units
if [ "$repeat_length_remainder" -ne 0 ]; then
    >&2 echo "ERROR: contig length is not a multiple of the number of repeat units"
    exit 1;
fi

# check that window_force_boundary is within repeat_length if window_force_boundary is not NA
if [ ! "$window_force_boundary" == "NA" ] && [ ! "$window_force_boundary" -lt "$repeat_length" ]; then
    >&2 echo "ERROR: window_force_boundary is not within repeat_length"
    exit 1;
fi

# extract the origin coordinates from the bed file bed_file_with_features
# =======================================================================

if [ ! -f "$bed_file_with_features" ]; then
    >&2 echo "ERROR: bed file with features does not exist"
    exit 1;
fi

temp_bed_file_with_features=$(mktemp)
< "$bed_file_with_features" grep -E -v '^browser|^track|^#' > "$temp_bed_file_with_features"

while IFS=$'\t' read -r _ col2 col3 col4 _ _; do
  if [[ "$col4" == "ARS1200-1" ]]; then
    origin_start="$col2"
    origin_end="$col3"
    break
  fi
done < "$temp_bed_file_with_features"

# ensure that window_force_boundary is not within origin_start to origin_end
if [ ! "$window_force_boundary" == "NA" ] && [ "$window_force_boundary" -ge "$origin_start" ] &&\
       [ "$window_force_boundary" -le "$origin_end" ]; then
  >&2 echo "ERROR: window_force_boundary is within origin_start to origin_end"
  exit 1;
fi

# make temporary bed files along plus, minus, and "all" strands
# =============================================================

bed_file_plus=$(mktemp)
bed_file_minus=$(mktemp)
bed_file_all=$(mktemp)

# send a bed file entry of one rDNA unit to these files
{
if [ "$window_force_boundary" == "NA" ]; then
  echo -e "$contig_name\t0\t$repeat_length\tFULL_CONTIG_PLUS\t0\t+\t$((repeat_length - 1))"
else
  echo -e "$contig_name\t0\t$repeat_length\tFULL_CONTIG_PLUS\t0\t+\t$window_force_boundary"
fi
} > "$bed_file_plus"

# flip the strand of the bed file to create the minus strand bed file
# shellcheck disable=SC1010
# shellcheck disable=SC2016
mlr --itsv --otsv --implicit-csv-header --headerless-csv-output put 'if($6=="+"){$6="-"}elif($6=="-"){$6="+"}' then cat\
  "$bed_file_plus" > "$bed_file_minus"

# Prepare bed file along reference for any data that cannot be mapped to a strand
< "$bed_file_plus" awk -v OFS='\t' '{print $1, $2, $3, $4, $5, ".", $7}' > "$bed_file_all"

# make windows on the plus, minus, all bed files on the reference, and count bases per window
# ===========================================================================================

bed_files=("$bed_file_plus" "$bed_file_minus" "$bed_file_all")
output_bed_files=("reference_windows_plus.bed" "reference_windows_minus.bed" "reference_windows_all.bed")

# Initialize a counter
counter=0

# declare bases and their respective columns in our counting script
declare -a bases=("adenosine" "cytidine" "guanosine" "thymidine" "all")
declare -a columns=("7" "8" "9" "10" "11")

# make windows on the plus, minus, all bed files
for bed_file in "${bed_files[@]}"; do

  # make windows
  # IMPORTANT NOTE: if you are going to make windows in thymidines instead of nucleotides in general,
  # you've to think very carefully about what windows on the reference mean.
  python make_windows_given_fasta_and_bed.py -i "$bed_file" -g "$ref_fasta"\
    -n "$window_number_bases" -b "$window_base" -d FB |\
    sort -k1,1 -k2,2n | tee "$output_dir"/"${output_bed_files[$counter]}"_temp |\
      awk 'BEGIN {OFS="\t"} {print $1, $2, $3}' | insert_calling_script_header "$@" \
        > "$output_dir"/"${output_bed_files[$counter]}"

  # count number of all types of bases along windows
  python count_bases_windows_given_fasta_and_bed.py -i "$output_dir"/"${output_bed_files[$counter]}"_temp -g "$ref_fasta" >\
    "$output_dir"/"${output_bed_files[$counter]}"_temp_count

  # separate information by base
  for ((i = 0; i < ${#bases[@]}; i++)); do

    base="${bases[$i]}"
    column="${columns[$i]}"

    < "$output_dir"/"${output_bed_files[$counter]}"_temp_count \
      awk -v col="$column" 'BEGIN{OFS=" "}{print $1, $2, $3, $(col)}' |\
        insert_calling_script_header "$@" > "$output_dir"/"${output_bed_files[$counter]}"_"$base".bedgraph

  done

  # remove temporary file
  rm "$output_dir"/"${output_bed_files[$counter]}"_temp
  rm "$output_dir"/"${output_bed_files[$counter]}"_temp_count

  # Increment the counter
  counter=$((counter+1))

done

# convert pause file and separate pauses along + and - strands
# ============================================================

< "$pause_file" python convert_pause_file_to_bed.py --discardUnkeptPauses |\
  sed "1i# from $pause_file" | insert_calling_script_header "$@" \
    > "$output_dir"/pauseFile_all.bed

< "$pause_file" python convert_pause_file_to_bed.py --discardUnkeptPauses |\
  awk '{if ($6 == "+") {print $0}}' | sed "1i# from $pause_file" | insert_calling_script_header "$@" \
    > "$output_dir"/pauseFile_plus.bed

< "$pause_file" python convert_pause_file_to_bed.py --discardUnkeptPauses |\
  awk '{if ($6 == "-") {print $0}}' | sed "1i# from $pause_file" | insert_calling_script_header "$@" \
    > "$output_dir"/pauseFile_minus.bed

# convert pause file to lead/lag bed file
< "$pause_file" python convert_pause_file_to_bed.py --discardUnkeptPauses --LeadLagToPlusMinus |\
  awk '{if ($6 == "+") {print $0}}' | sed "1i# from $pause_file" | insert_calling_script_header "$@" \
    > "$output_dir"/pauseFile_lead.bed

< "$pause_file" python convert_pause_file_to_bed.py --discardUnkeptPauses --LeadLagToPlusMinus |\
  awk '{if ($6 == "-") {print $0}}' | sed "1i# from $pause_file" | insert_calling_script_header "$@" \
    > "$output_dir"/pauseFile_lag.bed

# convert pause file to L/R bed file
< "$pause_file" python convert_pause_file_to_bed.py --discardUnkeptPauses --LRtoPlusMinus |\
  awk '{if ($6 == "-") {print $0}}' | sed "1i# from $pause_file" | insert_calling_script_header "$@" \
    > "$output_dir"/pauseFile_left.bed

< "$pause_file" python convert_pause_file_to_bed.py --discardUnkeptPauses --LRtoPlusMinus |\
  awk '{if ($6 == "+") {print $0}}' | sed "1i# from $pause_file" | insert_calling_script_header "$@" \
    > "$output_dir"/pauseFile_right.bed

# take care of periodic coordinates for pause files and nascent and parent bed files, and get coverage
# ====================================================================================================

# create an array of input and output files
input_bed_files=()
output_bed_files=()
reference_files=()
output_bedgraphs=()

# pause intervals are much smaller than window sizes as pause intervals are 1 bp long.
# alignment-coordinate intervals are much larger than window sizes as these are several kb long.
# so, in one case, we need to count coverage of the feature by the window and in the other case the coverage
# of the window by the feature. we accomplish this through an inversion flag (-i) in the coverage script.
inversion_flag=()

# populate the input output arrays for different pause bed files and parent and nascent bed files if available
for suffix in plus minus lead lag left right all; do
  input_bed_files+=("$output_dir"/pauseFile_"$suffix".bed)
  output_bed_files+=("$output_dir"/pauseFile_"$suffix"_fix_periodic.bed)
  output_bedgraphs+=("$output_dir"/pauseFile_"$suffix"_fix_periodic.bedgraph)
  if [ "$suffix" == "plus" ] || [ "$suffix" == "minus" ]; then
    reference_files+=("$output_dir"/reference_windows_"$suffix".bed)
  else
    reference_files+=("$output_dir"/reference_windows_all.bed)
  fi
  inversion_flag+=("-i")
done

if [ -f "$nascent_bed" ]; then
  input_bed_files+=("$nascent_bed")
  output_bed_files+=("$output_dir"/nascent_fix_periodic.bed)
  output_bedgraphs+=("$output_dir"/nascent_fix_periodic.bedgraph)
  reference_files+=("$output_dir"/reference_windows_all.bed)
  inversion_flag+=("")
fi

if [ -f "$parent_bed" ]; then
  input_bed_files+=("$parent_bed")
  output_bed_files+=("$output_dir"/parent_fix_periodic.bed)
  output_bedgraphs+=("$output_dir"/parent_fix_periodic.bedgraph)
  reference_files+=("$output_dir"/reference_windows_all.bed)
  inversion_flag+=("")
fi

# account for periodic coordinates, and get bedgraphs of coverage
for counter in "${!input_bed_files[@]}"; do

  bash convert_normal_bed_coords_to_periodic_bed_coords_rDNA.sh "${input_bed_files[$counter]}" \
    "$ref_fasta_index" "$contig_name" "$num_repeat_units" intensive > "${output_bed_files[$counter]}"

  # shellcheck disable=SC2086
  bash convert_bed_to_coverage_given_windows.sh ${inversion_flag[$counter]} "${reference_files[$counter]}" \
    "${output_bed_files[$counter]}" "$ref_fasta_index" | insert_calling_script_header "$@" \
      > "${output_bedgraphs[$counter]}"

done


# normalize coverage of various pause bed files by coverage of reads if nascent and parental bed files are available
# ==================================================================================================================

if [ -f "$nascent_bed" ] && [ -f "$parent_bed" ]; then

  # add parent and nascent bedgraphs
  python perform_binary_bedgraph_operation.py\
    "$output_dir"/nascent_fix_periodic.bedgraph\
    "$output_dir"/parent_fix_periodic.bedgraph add > "$output_dir"/all_fix_periodic.bedgraph

  # normalize pause bedgraphs by all bedgraph count and further by base count
  for suffix in left right lead lag plus minus all; do
    python divide_bedgraph_values.py "$output_dir"/pauseFile_"$suffix"_fix_periodic.bedgraph \
      "$output_dir"/all_fix_periodic.bedgraph |\
        insert_calling_script_header "$@" |\
        sed '1i# values multiplied by 1 million for sake of IGV view' |\
          awk '/^#/{print} !/^#/ {$4 = $4*1000000; print}' > \
            "$output_dir"/pauseFile_"$suffix"_fix_periodic_norm_1milX.bedgraph

    python divide_bedgraph_values.py "$output_dir"/pauseFile_"$suffix"_fix_periodic_norm_1milX.bedgraph \
      "${output_dir}/reference_windows_all.bed_all.bedgraph" |\
        insert_calling_script_header "$@" |\
        sed '1i# values multiplied by 1 million for sake of IGV view' >\
            "$output_dir"/pauseFile_"$suffix"_fix_periodic_norm_1milX_norm_baseCount.bedgraph
  done

fi

# remove temporary files
# ======================

rm "$bed_file_plus" "$bed_file_minus" "$bed_file_all" "$temp_bed_file_with_features"