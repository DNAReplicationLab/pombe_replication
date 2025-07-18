#!/bin/bash

#SBATCH --mem-per-cpu=10G
#SBATCH -c 10
#SBATCH -p ei-medium
#SBATCH -J getMnBPAcrossBedAfterSplit
#SBATCH --mail-type=END,FAIL
#SBATCH --time 23:59:59
#SBATCH --constraint=""

# goal
# ----

# Tile each entry in a bed file, putting all first tiles in one file, all second tiles in another, etc.
# Then, output mean BrdU per tile file.
# Use the script if for instance you want to collapse all pauses together and get mean BrdU vs coordinate
# along collapsed pauses.

# usage
# -----
# sbatch <program_name.sh> $mod_bam $bed_file $window_size [$tmp_dir] [$align_by] [$bed_filter_options] [$bed_op_dir] \
#         [minus_start_bp] [plus_end_bp]
# NOTE: [] means optional argument. If you want to specify an optional argument, then you have to specify all
#       optional arguments to the left of it. You must also remove the square brackets.
# NOTE: \ means the command continues in the next line.
# mod_bam: .mod.bam file containing analogue modification probabilities
# bed_file: bed file containing reads to tile, must have read_id in the fourth column i.e. as the name.
#           All intervals must be of the same size and the size must be divisible by window_size and by 2.
#           If there are multiple intervals with the same read id, then they must not overlap.
#           NOTE: If you made a bed file by growing some intervals, and the intervals get clipped because they
#                 extend beyond the chromosome, then you must discard those intervals before running this script as they
#                 will cause an error due to the uneven size of the intervals.
# window_size: window size of the tile in bp. WARNING: if you expect a large number of tiles, then adjust number
#              of cores requested accordingly.
# tmp_dir: (optional, defaults to /tmp) temporary directory to store intermediate files
# align_by: (optional, defaults to head-to-tail) before tiling, align bed entries head-to-tail or tail-to-head or
#           by increasing-ref or by decreasing-ref. See the section below for more information.
# bed_filter_options: (optional) options to filter bed file, default "forks". If the bed file came from conversion
#                     from a pause file in our format, it would contain information about fork calls and alignment
#                     coordinates. So we can use that to reject bed intervals that lie outside the fork or the
#                     alignment. So, please set this to "forks" (default) or "alignments". If the bed file was
#                     generated from a different source, then this option has no effect. We basically search for
#                     a column in the bed file that is not part of the standard bed format and use that to filter
#                     out intervals. If that column does not exist, then we do not filter out any intervals.
# bed_op_dir: (optional) directory to store output bed files. If not specified, then the output bed files are
#             discarded. Please note that these are files output in addition to the table that is written to
#             standard output, which is unaffected by this option.
# minus_start_bp: (optional) subtract so many bp from the start of each interval before tiling. default 0.
#                 user must ensure that the final interval size is even and is a multiple of window_size.
# plus_end_bp: (optional) add so many bp to the end of each interval before tiling. default 0.
#              user must ensure that the final interval size is even and is a multiple of window_size.

# Many columns are output.

# Meaning of align_by
# -------------------
# Every bed entry has a direction associated with i.e. + or -, which means 5' to 3' or 3' to 5' respectively.
# So for tiling, we have four choices, we can align bed entries head-to-tail, tail-to-head,
# by increasing-ref or by decreasing-ref.
# In the above options, the "head" means the head of the arrow and the "tail" means the tail of the arrow, where
# the arrow is the direction of the bed entry, meaning it points from 5' to 3' (+) or 3' to 5' (-).
# An advanced usage is if you want to align by fork direction, you create the bed file with the fork direction
# represented in the sixth column as + or -. Then, you can use head-to-tail or tail-to-head to align by fork
# direction. You can create such a bed file from our pause file using the script convert_pause_file_to_bed.py
# using suitable options.

# Example for tiling
# ------------------
# Let's say the entries are (for ease of readability, I am using read1 and read2 as the name of the reads, instead
# of UUIDs):
# chr1 100 200 read1 + 0
# chr1 300 400 read2 - 0
# if you choose a window size of 10 and tiling of head-to-tail, then the contents of tile 1 are:
# chr1 190 200 read1 + 0
# chr1 300 310 read2 - 0
# and the contents of tile 2 are:
# chr1 180 190 read1 + 0
# chr1 310 320 read2 - 0
# and so on.
# this is because the "head" of read1 is at 200 and the "head" of read2 is at 300.
# if you chose a tiling of increasing-ref, then the contents of tile 1 are:
# chr1 100 110 read1 + 0
# chr1 300 310 read2 - 0
# and the contents of tile 2 are:
# chr1 110 120 read1 + 0
# chr1 310 320 read2 - 0
# and so on.

# stop execution if any command fails
set -e

# need at least three arguments
if [ $# -lt 3 ]
then
    >&2 echo "Goal: tile each bed entry and get mean brdu per tile."
    >&2 echo "      for more information, please see the header of the script."
    >&2 echo "usage: sbatch <program_name.sh> <mod_bam> <bed_file> <window_size> [tmp_dir] [align_by] \ "
    >&2 echo "         [bed_filter_options] [bed_op_dir] [minus_start_bp] [plus_end_bp]"
    >&2 echo "mod_bam: .mod.bam file containing analogue modification probabilities"
    >&2 echo "bed_file: bed file containing reads to tile, must have read_id in the fourth column i.e. as the name."
    >&2 echo "          Please read the header of the script for more information about the bed file format."
    >&2 echo "window_size: window size of the tile in bp. WARNING: if you expect a large number of tiles, "
    >&2 echo "              then adjust number of cores requested in the script accordingly."
    >&2 echo "tmp_dir: (optional, defaults to /tmp) temporary directory to store intermediate files"
    >&2 echo "align_by: (optional, defaults to head-to-tail) before tiling, align bed entries head-to-tail or "
    >&2 echo "          tail-to-head or by increasing-ref or by decreasing-ref."
    >&2 echo "bed_filter_options: (optional) options to filter bed file, default \"forks\". "
    >&2 echo "                    Please read the header of the script for more information."
    >&2 echo "bed_op_dir: (optional) directory to store output bed files. "
    >&2 echo "            Please read the header of the script for more information."
    >&2 echo "minus_start_bp: (optional) subtract so many bp from the start of each interval before tiling. default 0."
    >&2 echo "                user must ensure that the final interval size is even and is a multiple of window_size."
    >&2 echo "plus_end_bp: (optional) add so many bp to the end of each interval before tiling. default 0."
    >&2 echo "             user must ensure that the final interval size is even and is a multiple of window_size."
    exit 1
fi

# load packages
source load_package.sh -python -samtools -jq -miller -bedtools

# load configuration
source config.sh

# accept input arguments
mod_bam=$1
bed_file=$2
window_size=$3
temp_dir_root=${4:-"${config[scratchDir]:-}"/tmp}
align_by=${5:-head-to-tail}
bed_filter_options=${6:-forks}
bed_op_dir=${7:-}
minus_start_bp=${8:-0}
plus_end_bp=${9:-0}

# if bed_filter_options is not one of forks or alignments, then set it to forks
if [ "$bed_filter_options" != "forks" ] && [ "$bed_filter_options" != "alignments" ]; then
    bed_filter_options=forks
fi

# use full paths for mod_bam and bed_file
mod_bam=$(readlink -f "$mod_bam")
bed_file=$(readlink -f "$bed_file")

# check if temp_dir_root exists
if [ ! -d "$temp_dir_root" ]; then
    >&2 echo "Error: $temp_dir_root does not exist."
    exit 1
fi

# create a temporary directory relative to temp_dir_root
temp_dir=$(mktemp -d -p "$temp_dir_root")

# check that mod_bam exists
if [ ! -f "$mod_bam" ]; then
    >&2 echo "Error: $mod_bam does not exist."
    exit 1
fi

# check that the mod bam index exists
if [ ! -f "$mod_bam.bai" ]; then
    >&2 echo "Error: $mod_bam.bai does not exist."
    exit 1
fi

# check that bed_file exists
if [ ! -f "$bed_file" ]; then
    >&2 echo "Error: $bed_file does not exist."
    exit 1
fi

# validate bed file
if [ ! "$(< "$bed_file" python validate_bed_format.py --require-uuid --allow-float-score)" == "valid" ]; then
    >&2 echo "Error: $bed_file is not in the correct format."
    exit 1;
fi

# check that window_size is a positive integer
if ! [[ "$window_size" =~ ^[0-9]+$ ]]; then
    >&2 echo "Error: $window_size is not a positive integer."
    exit 1
fi

# check that window_size is greater than 0
if [ "$window_size" -le 0 ]; then
    >&2 echo "Error: $window_size is not greater than 0."
    exit 1
fi

# check that the minus_start_bp and plus end bp are greater than or equal to 0
for bp in "$minus_start_bp" "$plus_end_bp"; do
  if ! [[ "$bp" =~ ^[0-9]+$ ]]; then
      >&2 echo "Error: $bp is not a positive integer."
      exit 1
  fi
  if [ "$bp" -lt 0 ]; then
      >&2 echo "Error: $bp is not greater than or equal to 0."
      exit 1
  fi
done

# check that align_by is one of head-to-tail, tail-to-head, increasing-ref, decreasing-ref
if [ "$align_by" != "head-to-tail" ] && [ "$align_by" != "tail-to-head" ] \
  && [ "$align_by" != "increasing-ref" ] && [ "$align_by" != "decreasing-ref" ]; then
    >&2 echo "Error: $align_by is not one of head-to-tail, tail-to-head, increasing-ref, decreasing-ref."
    exit 1
fi

# check if bed file has self-intersections with same names
# ========================================================
# first, extract first six columns from the bed file
tmp_bed_file=$(mktemp "$temp_dir"/XXXXXXXXXX)
< "$bed_file" grep -E -v '^browser|^track|^#' |\
  awk 'BEGIN{OFS="\t"}{if($2!=$3){print $1,$2,$3,$4,$5,$6}}' > "$tmp_bed_file"
# shellcheck disable=SC1010
# shellcheck disable=SC2016
number_of_intersections=\
$(bedtools intersect -a "$tmp_bed_file" -b "$tmp_bed_file" -wo |\
   mlr --tsv --skip-comments --implicit-csv-header --headerless-csv-output \
   filter '($4 == $10)' then count)
number_of_lines=$(wc -l < "$tmp_bed_file")
if [ "$number_of_intersections" -gt "$number_of_lines" ]; then
  >&2 echo "Error: $bed_file has self-intersections with same names."
  rm "$tmp_bed_file"
  exit 1
fi

# check that all intervals are of the same size and the size is divisible by window_size and by 2
# ===============================================================================================
# check that all intervals are of the same size
if [ "$(awk '{print $3-$2}' "$tmp_bed_file" | sort -u | wc -l)" -ne 1 ]; then
  >&2 echo "Error: all intervals in $bed_file must be of the same size."
  rm "$tmp_bed_file"
  exit 1
fi

# get the size of the first interval and grow it if needed
interval_size=$(awk '{print $3-$2}' "$tmp_bed_file" | head -n 1)
interval_size_grown=$((interval_size + minus_start_bp + plus_end_bp))

# check that the size is divisible by window_size and by 2 and is greater than window_size
if [ "$((interval_size_grown % window_size))" -ne 0 ] || [ "$((interval_size_grown % 2))" -ne 0 ] || \
[ "$interval_size_grown" -le "$window_size" ]; then
  >&2 echo "Error: the size of all intervals in $bed_file must be divisible by window_size and by 2 "
  >&2 echo "       and greater than window_size."
  rm "$tmp_bed_file"
  exit 1
fi

# subset mod bam file by reads in bed file
# ========================================
# extract 4th column of bed file
grep -v '^#' "$bed_file" | awk '{print $4}'  > "$temp_dir/bed_4th_col.txt"

# use samtools to extract reads from mod_bam that are in bed file
samtools view -@ 15 -h -N "$temp_dir/bed_4th_col.txt" "$mod_bam" > "$temp_dir/mod_bam_subset.bam"

# sort and index the subset mod bam file
samtools sort -@ 15 -o "$temp_dir/mod_bam_subset_sorted.bam" "$temp_dir/mod_bam_subset.bam"
samtools index "$temp_dir/mod_bam_subset_sorted.bam"

# split bed file
# ===============
# NOTE: we add a ridiculously large number here because we do not want negative coordinates.
# We will subtract this number later.
< "$bed_file" grep -v '^#' |\
  awk -v a="$minus_start_bp" -v b="$plus_end_bp" -v OFS="\t" '{$2=$2-a+10000000000; $3=$3+b+10000000000; print $0}' |\
    python split_each_bed_file_entry_along_length.py --prefix "$temp_dir/split" \
      --num-bases "$window_size" --align-by "$align_by" > /dev/null

# get mean BrdU per tile and store it
# ====================================
thread_count=10
current_thread=0
for file in "$temp_dir/split_"*".bed"; do

  # file with BrdU signal
  bed_with_brdu_signal=${file%.bed}.mod.bed

  # reject any interval that lies outside the corresponding fork and get statistics.
  # we subtract the large number as negative coordinates will get naturally trimmed by our rejection criterion above.
  # shellcheck disable=SC2094
  < "$file" awk -v OFS="\t" '{$2=$2-10000000000; $3=$3-10000000000; print $0}' |\
  python filter_pause_file_bed_intervals_outside.py --"$bed_filter_options" |\
    bash get_mean_data_from_modBAM_using_bed.sh - "$temp_dir/mod_bam_subset_sorted.bam" | tee "$bed_with_brdu_signal" |\
      mlr --itsv --ojson --skip-comments --implicit-csv-header stats1 -a \
        mean,stddev,count,min,max,p5,p10,p15,p20,p25,p30,p35,p40,p45,p50,p55,p60,p65,p70,p75,p80,p85,p90,p95 -f 5 |\
          jq -r '.[]' |\
            sed '2i\"input_bed_file\": \"'"$(basename "$file")"'\",' |\
              mlr --ijson --otsv --headerless-csv-output cat > "$file".stats &

  # increment thread count
  current_thread=$((current_thread + 1))

  # check if we have reached the maximum number of threads
  if [ "$current_thread" -eq "$thread_count" ]; then
    # wait for all jobs to finish
    wait;
    # reset thread count
    current_thread=0
  fi

done

# wait for all jobs to finish
wait;

# print results
# =============
echo -e "file\tmean\tsd\tcount\tmin\tmax\tp5\tp10\tp15\tp20\tp25\tp30\tp35\tp40\tp45\tp50\tp55\tp60\tp65\tp70\tp75\tp80\tp85\tp90\tp95"

for file in "$temp_dir/split_"*".bed.stats"; do
  cat "$file"
done

# store bed files if requested
# =============================
# move bed files
if [ -n "$bed_op_dir" ]; then
  mkdir -p "$bed_op_dir"
  mv "$temp_dir/split_"*".mod.bed" "$bed_op_dir"
fi

# bed files are strictly speaking not in the correct format as a strand column is missing.
# so, edit the bed files and add a strand column '.' to them
for file in "$bed_op_dir/split_"*".mod.bed"; do
  < "$file" grep -E -v '^browser|^track|^#' | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,"."}' > "$file.tmp"
  mv "$file.tmp" "$file"
done

# concatenate all bed files into one
find "$bed_op_dir" -name 'split_*.mod.bed' -type f  | sort -V | xargs cat > "$bed_op_dir/all_tiles.mod.bed"

# clean up
# ========
# remove temporary directory
rm -r "$temp_dir"