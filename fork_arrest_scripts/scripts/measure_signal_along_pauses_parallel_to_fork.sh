#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 2
#SBATCH -p ei-medium
#SBATCH -J getSigAlongPauses
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# goal
# -----
# Aggregate some signal available in a bedgraph in the neighbourhood of detected pauses and average it.

# An example
# ----------
# It is easier to understand this script with an example.
# Let's say we have the density of some component like a nucleosome in a bedgraph.
# We want to look in a 100 bp neighbourhood of each detected pause, and plot the nucleosome density
# from x = -50 to x = +50 averaged over all pauses in bins of size 5 bp.
# Moreover, we want the direction of increasing x to be the fork direction i.e. flip neighbourhoods depending
# on whether pauses are on left-moving or right-moving forks.
# Then, this script will do that for you.

# usage
#------
# bash measure_signal_along_pauses_parallel_to_fork.sh $pause_file $bedgraph_file $fasta_fai $output_dir $slop_size\
#   $bin_size [$filter_option] [$fasta_file] [$pause_sens_bedgraph_prefix]
# NOTE: [] means optional argument. If you want to specify one, remove the [] and write the argument within.
#       If you don't want to specify all of them, just ignore anything beginning with [].
#       If you want to specify some of them, you must specify the ones before it.
# $pause_file: pauses in our tab-separated format with column names.
#              * Basically, we need at least two columns: detectIndex, pauseSite.
#              * Any column that starts with keep is set to either True or False, and a pause is considered valid
#                only if all keep columns are True or if no keep columns are present in the file at all.
#              * detectIndex has the format readID_contig_start_end_strand_direction_startFork_endFork.
#                Names are self-explanatory. strand=fwd/rev; direction=L/R; start < end; startFork < endFork
#                irrespective of strand and fork direction.
# $bedgraph_file: bedgraph file with signal to be measured.
#                 NOTE: If you want missing intervals to be counted as zero, you must explicitly add them to the
#                       bedgraph file with signal 0, otherwise missing intervals will be treated as missing.
# $fasta_fai: fasta index file. Comes with the fasta file or can be generated using samtools faidx.
# $output_dir: directory where the output will be written.
# $slop_size: number of base pairs to the left and right of the pause to measure the signal.
#             In our example above, this is 50.
#             NOTE: pauses are 1 bp wide, so given a slop_size you'll end up with a neighbourhood of 2*slop_size+1 bp.
# $bin_size: number of base pairs to bin the signal.
#            In our example above, this is 5.
#            NOTE: For a division with equal bases per bin, 2*slop_size+1 must be divisible by bin_size.
# $filter_option: Option to filter forks for the pauses. left/right/fwd/rev/lead/lag/all are the options.
#                 Defaults to all. Left/right means left/right-moving forks. Fwd/rev means fwd/rev alignments.
#                 Lead/lag means forks on the leading/lagging strand. All means all forks.
# $fasta_file: (optional, default unused) fasta file. Used to plot AT count in the neighbourhood of pauses.
#              We could have inferred this from the fasta index file, but historically this script used only the
#              fasta index file and produced no such plot, so we want to keep the behaviour consistent.
#              If not specified, the AT count plot will not be produced.
#              This is a control to see if the feature we are measuring is real or not (see the description of the
#              parameter below).
# $pause_sens_bedgraph_prefix: (optional, default unused) prefix of files containing our null hypothesis for pause
#                              count, produced by the pipeline.
#                              Used to plot pause sensitivity in the neighbourhood of pauses.
#                              If not specified, this plot is not produced.
#                              We need controls for the calculations here.
#                              This and the AT count plot are two such controls.
#                              We are asking for a prefix here, as there are several sensitivity tracks called 'all',
#                              'left', 'left.minus' etc. and we need to be able to use the correct file depending on the
#                              filter_option. So, if the prefix is say /a/b/c/pause_sensitivity, then the script will
#                              look for /a/b/c/pause_sensitivity.all.bedgraph, /a/b/c/pause_sensitivity.left.bedgraph
#                              etc. as needed.

# outputs
# -------
# A bunch of bed files, a tab-separated file, and a plot are sent to the output directory.
# We are not going to describe them in detail here. Consult the respective scripts below for more details.
# The overall details are: We split intervals corresponding to pauses by length along the fork direction
#                          to form the bed files. Then, we measure a mean and sd per bed file and send it to the
#                          tab-separated file, and plot it.

# stop execution if any command fails
set -e

# load packages
source load_package.sh -python -bedtools -seqtk

# load configuration
source config.sh

# set temporary directory
tmpDir="${config[scratchDir]:-}"/tmp/temp$(openssl rand -hex 6)
mkdir -p "$tmpDir"

# check that the correct number of arguments were provided
if [ "$#" -lt 6 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash measure_signal_along_pauses_parallel_to_fork.sh \$pause_file \$bedgraph_file \$fasta_fai \$output_dir \$slop_size \$bin_size \$filter_option"
    >&2 echo "\$pause_file: pauses in our tab-separated format."
    >&2 echo "\$bedgraph_file: bedgraph file with signal to be measured."
    >&2 echo "\$fasta_fai: fasta index file. Comes with the fasta file or can be generated using samtools faidx."
    >&2 echo "\$output_dir: directory where the output will be written."
    >&2 echo "\$slop_size: number of base pairs to the left and right of the pause to measure the signal."
    >&2 echo "       NOTE: pauses are 1 bp wide, so given a slop_size you'll end up with a neighbourhood of 2*slop_size+1 bp."
    >&2 echo "\$bin_size: number of base pairs to bin the signal."
    >&2 echo "       NOTE: For a division with equal bases per bin, 2*slop_size+1 must be divisible by bin_size."
    >&2 echo "\$filter_option: Option to filter forks for the pauses. left/right/fwd/rev/lead/lag/all are the options. Defaults to all."
    >&2 echo "\$fasta_file: (optional, default unused) fasta file. Used to plot AT count in the neighbourhood of pauses."
    >&2 echo "              If not specified, the AT count plot will not be produced."
    >&2 echo "\$pause_sens_bedgraph_prefix: (optional, default unused) our null hypothesis for pause count, produced by"
    >&2 echo "                              the pipeline. Used to plot pause sensitivity in the neighbourhood of pauses."
    >&2 echo "It is highly recommended to read the comments in the script as the description printed here is terse."
    exit 1;
fi

# parse arguments and check they are valid
# ----------------------------------------
echo "INFO: parsing arguments and checking if files exist/arguments are valid etc."
pause_file="$1"
bedgraph_file="$2"
fasta_fai="$3"
output_dir="$4"
slop_size="$5"
bin_size="$6"
filter_option=${7:-"all"}
fasta_file=${8:-/dev/null}
pause_sens_bedgraph_prefix=${9:-/dev/null}

# check if the files exist
if [ ! -f "$pause_file" ]; then
    >&2 echo "ERROR: pause file $pause_file does not exist"
    exit 1;
fi

if [ ! -f "$bedgraph_file" ]; then
    >&2 echo "ERROR: bedgraph file $bedgraph_file does not exist"
    exit 1;
fi

if [ ! -f "$fasta_fai" ]; then
    >&2 echo "ERROR: fasta index file $fasta_fai does not exist"
    exit 1;
fi

# check if the output directory exists or else create it
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

output_dir=$(realpath "$output_dir")

# check that slop size and bin size are positive integers
if ! [[ "$slop_size" =~ ^[0-9]+$ ]]; then
    >&2 echo "ERROR: slop size $slop_size is not a positive integer"
    exit 1;
fi

if ! [[ "$bin_size" =~ ^[0-9]+$ ]]; then
    >&2 echo "ERROR: bin size $bin_size is not a positive integer"
    exit 1;
fi

if [ "$slop_size" -lt 1 ]; then
    >&2 echo "ERROR: slop size $slop_size is not a positive integer"
    exit 1;
fi

if [ "$bin_size" -lt 1 ]; then
    >&2 echo "ERROR: bin size $bin_size is not a positive integer"
    exit 1;
fi

# check that pause file and bedgraph file are valid
# -------------------------------------------------
# check that the bedgraph file, pause file, and fasta index files have multiple lines
if [ "$(< "$pause_file" grep -c -v "^#")" -lt 2 ]; then
    >&2 echo "ERROR: pause file $pause_file has less than 2 lines"
    exit 1;
fi
if [ "$(< "$bedgraph_file" grep -c -E -v "^#|^browser|^track")" -lt 2 ]; then
    >&2 echo "ERROR: bedgraph file $bedgraph_file has less than 2 lines"
    exit 1;
fi
if [ "$(< "$fasta_fai" grep -c -v "^#")" -lt 2 ]; then
    >&2 echo "ERROR: fasta index file $fasta_fai has less than 2 lines"
    exit 1;
fi

# validate pause file
if [ ! "$(< "$pause_file" python validate_pause_format.py)" == "valid" ]; then
  >&2 echo "ERROR: pause_file $pause_file is not valid"
  exit 1;
fi

# we are converting the bedgraph to a bed file and checking if it is valid
bash convert_bedgraph_to_bed.sh -u "$bedgraph_file" -f "$fasta_fai"  > "$tmpDir"/bedgraph.bed

if [ ! "$(< "$tmpDir"/bedgraph.bed python validate_bed_format.py --never-zero-base)" == "valid" ]; then
    >&2 echo "Error: something wrong with bedgraph file $bedgraph_file"
    exit 1;
fi

# check that the bed file is valid
if [ ! "$(< "$tmpDir"/bedgraph.bed python validate_bed_against_fai.py "$fasta_fai" )" == "valid"  ]; then
  >&2 echo "Error: something wrong with bedgraph file $bedgraph_file"
  exit 1;
fi

rm "$tmpDir"/bedgraph.bed

# convert pause file to bed (using fork direction for strand sign) and grow pauses by the given size
# --------------------------------------------------------------------------------------------------
echo "INFO: converting pause file to bed and slopping it by $slop_size bp and filtering by $filter_option"
echo "INFO: please note that an invalid prefix results in a zero-sized output and could cause downstream errors"
slopped_pause_file="$tmpDir"/slopped_pause_file.bed
# FLAG: LEAD LAG DISTINCTION
< "$pause_file" python convert_pause_file_to_bed.py --LRtoPlusMinus --discardUnkeptPauses \
    --outputBed6Plus1WherePlus1isAlignStrand |\
  bedtools slop -i - -g "$fasta_fai" -b "$slop_size" |\
awk -F '\t' -v filter="$filter_option" '
{
    if (filter == "left") { if ($6 == "-") print }
    else if (filter == "right") { if ($6 == "+") print }
    else if (filter == "fwd") { if ($7 == "+") print }
    else if (filter == "rev") { if ($7 == "-") print }
    else if (filter == "lead") { if (($6 == "+" && $7 == "+") || ($6 == "-" && $7 == "-")) print }
    else if (filter == "lag") { if (($6 == "+" && $7 == "-") || ($6 == "-" && $7 == "+")) print }
    else if (filter == "all") { print }
}' | sed '1i# -/+ in sixth column means left/right fork. Seventh column is the alignment strand' > "$slopped_pause_file"

# if there are no valid pauses, then exit
if [ "$(wc -l < "$slopped_pause_file")" -lt 1 ]; then
    >&2 echo "ERROR: no valid pauses found after filtering by $filter_option in $pause_file"
    exit 1;
fi

# split the slopped pause file by length along the fork direction
# ---------------------------------------------------------------
echo "INFO: splitting slopped pause file by length along the fork direction"

< "$slopped_pause_file" python split_each_bed_file_entry_along_length.py --prefix "$output_dir/split_pause_file" \
  --num-bases "$bin_size" --align-by "head-to-tail"

# Get and plot the signal
# ----------------------------
echo "INFO: getting and plotting the signal (in a background process)"
tmp_json_2="$tmpDir"/get_and_plot_mean_bedgraph_signal_per_bed_in_directory_options.json
{
  echo "{"
  echo "\"bed_sort_file_name_numerical\": true,"
  echo "\"output_file_prefix\": \"signal\"",
  echo "\"y_label_plot\": \"signal\""
  echo "}"
} > "$tmp_json_2"

bash get_and_plot_mean_bedgraph_signal_per_bed_in_directory.sh "$output_dir" "$bedgraph_file" "$output_dir"\
  "$tmp_json_2" "..i" &

# Get and plot the AT count if fasta file is provided
# ---------------------------------------------------
if [ -f "$fasta_file" ]; then
  echo "INFO: getting and plotting the AT count (in a background process)"
  tmp_json_3="$tmpDir"/get_and_plot_AT_count_per_bed_in_directory_options.json
  {
    echo "{"
    echo "\"bed_sort_file_name_numerical\": true,"
    echo "\"output_file_prefix\": \"AT_count\","
    echo "\"y_label_plot\": \"AT count\""
    echo "}"
  } > "$tmp_json_3"

  # compute the AT count in the neighbourhood of pauses
  tmp_bedgraph_file="$tmpDir"/$(openssl rand -hex 6)
  seqtk comp -r "$bedgraph_file" "$fasta_file" |\
    awk -v OFS=" " -v IFS="\t" '{print $1, $2, $3, $4 + $7}' > "$tmp_bedgraph_file"

  bash get_and_plot_mean_bedgraph_signal_per_bed_in_directory.sh "$output_dir" "$tmp_bedgraph_file" "$output_dir"\
    "$tmp_json_3" "..i" &
fi

# Get and plot the pause sensitivity if pause sensitivity bedgraph file is provided
# ---------------------------------------------------------------------------------
if [ -f "$pause_sens_bedgraph_prefix".all.bedgraph ]; then
  echo "INFO: getting and plotting the pause sensitivity (in a background process)"
  tmp_json_4="$tmpDir"/get_and_plot_mean_sens_per_bed_in_directory_options.json
  {
    echo "{"
    echo "\"bed_sort_file_name_numerical\": true,"
    echo "\"output_file_prefix\": \"pause_sensitivity\","
    echo "\"y_label_plot\": \"Pause sensitivity\""
    echo "}"
  } > "$tmp_json_4"

  # if filter is lead or lag, then we need to combine the bedgraph files
  # FLAG: LEAD LAG DISTINCTION
  tmp_bedgraph_file="$tmpDir"/$(openssl rand -hex 6)
  if [ "$filter_option" == "lead" ]; then
    echo "INFO: computing pause sensitivity for leading strand"
    file_1="$pause_sens_bedgraph_prefix".right.plus.bedgraph
    file_2="$pause_sens_bedgraph_prefix".left.minus.bedgraph
    python perform_binary_bedgraph_operation.py "$file_1" "$file_2" add > "$tmp_bedgraph_file"
  elif [ "$filter_option" == "lag" ]; then
    echo "INFO: computing pause sensitivity for lagging strand"
    file_1="$pause_sens_bedgraph_prefix".right.minus.bedgraph
    file_2="$pause_sens_bedgraph_prefix".left.plus.bedgraph
    python perform_binary_bedgraph_operation.py "$file_1" "$file_2" add > "$tmp_bedgraph_file"
  else
    cp "$pause_sens_bedgraph_prefix"."$filter_option".bedgraph "$tmp_bedgraph_file"
  fi

  # check that the bedgraph file has at least 2 lines
  if [ "$(< "$tmp_bedgraph_file" grep -c -E -v "^#|^browser|^track")" -lt 2 ]; then
      >&2 echo "ERROR: bedgraph file $tmp_bedgraph_file has less than 2 lines"
      exit 1;
  fi

  bash get_and_plot_mean_bedgraph_signal_per_bed_in_directory.sh "$output_dir" "$tmp_bedgraph_file" "$output_dir"\
    "$tmp_json_4" "..i" &
fi

# remove temporary directory after waiting for the background process to finish
# ------------------------------------------------------------------------------
echo "INFO: waiting for background processes to finish"
wait;
rm -rf "$tmpDir"