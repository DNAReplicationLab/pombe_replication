#!/bin/bash

# goal
# -----
# Count number of intervals of bed file A that occur in neighbourhood of intervals in bed file B.

# usage
#------
# cat $bed_file_B | bash count_bed_in_neighbourhood_of_bed.sh [-s|-S] [-d] $bed_file_A $fasta_fai $slop_bp $bin_size_bp \
#   $direction $op_dir [$A_score_filter]
# NOTE: \ is used above to indicate that the command continues on the next line.
# NOTE: [] indicates optional arguments. To specify them, remove the brackets and type the argument, in this case
#      -s or -S or -d or a number for $A_score_filter.
# NOTE: -s or -S is used to indicate whether three or six columns are used in the bed file A.
#       By default, we count un-stranded overlaps between A and B, but if -s or -S is used, we count stranded overlaps,
#       same or opposite strand respectively.
# NOTE: -d means do not merge the intervals of bed file A before counting. If file A contains genomic features that
#       overlap, then it may be useful to merge them before counting. But if file A is a bed file created from a pause
#       file, then we do want to count pauses that fall at the same bp as separate pauses.
# bed_file_A: bed file whose intervals are to be counted, see the section goal above.
#             6 columns are needed if -S or -s is used, otherwise 3 columns are enough.
#             We will perform a merge before counting, so overlapping or book-ended intervals are treated as one,
#             unless -d is used.
# bed_file_B: bed file in whose neighbourhood we will count the intervals of bed file A.
#             All intervals here must be 0 bp. 6 bed columns needed.
# fasta_fai: fasta index file for the reference genome.
# slop_bp: number of base pairs to extend the intervals of bed file B before counting the intervals of bed file A,
#          must be a multiple of bin_size_bp.
# bin_size_bp: size of the bins in base pairs to chop up the neighbourhood of B.
# direction: direction of the neighbourhood, can be 'head-to-tail', 'tail-to-head', 'increasing-ref', 'decreasing-ref'.
#            basically how to align intervals in B with each other to create the neighbourhood.
#            the 'head' and 'tail' above refer to arrow head and arrow tail where the arrow is the direction/strand
#            of the interval in the bed file.
# op_dir: output directory where files are sent, will be created if non-existent. Empty directories are preferred.
# A_score_filter: optional, if set, only intervals of bed file A with a score greater than or equal to this value
#                 are retained. Default is to ignore the score column whether it exists or not. To use this option,
#                 the A bed file must have at least six columns, the fifth of which is the score column.

# overall logic
# -------------
# Please also see details of logic below.
# 1. grow each interval in bed file B by +- slop bp.
# 2. divide the grown intervals into bins of size bin_size_bp: this gives us (2 * slop bp / bin_size_bp) bins.
#    put all intervals corresponding to each bin in its own bed file.
# 3. In each bed file above, count the number of intervals of bed file A that overlap with the bed file.
# 4. Output the counts in a tabular format.

# details of logic
# ----------------
# - for step 1, we want to grow intervals avoiding overlaps, so we want non-overlapping intervals (in a stranded sense)
#   in B all of zero length. So we can grow each interval, stopping as needed when we encounter intersections.
# - as a consequence of avoiding overlaps, each bed file generated in step two (which corresponds to each bin)
#   will have irregularly sized intervals in general, some intervals may even be 0 bp.
# - for step 3, for every bin, if there is any overlap between A and the bin, we will increment the count by 1.
#   We will do stranded counting if -s or -S is used, otherwise we will do unstranded counting.

# outputs
# -------
# - bed files corresponding to each bin are sent to the output directory.
# - an additional tabular file with counts per bin is also sent to the output directory. this file is tab-separated and
#   has four columns with column names (and comments starting with #): filename, distance, count, bp_in_bin.
#   filename is the name of the bin file, distance is the signed distance from the bin to the bed file B interval,
#   count is the number of intervals of bed file A that overlap with the bin, and bp_in_bin is the number of base pairs
#   in the bin.
# - the contents of this tabular file are also sent to stdout.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -bedtools -python -miller -jq

# load configuration
source config.sh

# set temporary directory
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

# complain if the number of arguments is less than 6
if [ $# -lt 6 ]; then
  >&2 echo "ERROR: Incorrect number of arguments"
  >&2 echo "Usage: cat bed_file_B | bash count_bed_in_neighbourhood_of_bed.sh [-s|-S] [-d] bed_file_A fasta_fai"\
" slop_bp bin_size_bp direction op_dir [A_score_filter]"
  >&2 echo "For meanings of arguments, see the section 'usage' in the script header."
  >&2 echo "The overall idea is that we grow file_B intervals by slop_bp, align them by the direction param, "
  >&2 echo "and count the number of intervals of file_A that overlap with each bin."
  exit 1;
fi

# parse options
# -------------
# set defaults
intersect_option=""
is_do_not_merge_A="false"
num_s_set=0

# parse arguments
while getopts ":sSd" opt; do
  case $opt in
    s)
      intersect_option="-s"
      num_s_set=$((num_s_set+1))
      ;;
    S)
      intersect_option="-S"
      num_s_set=$((num_s_set+1))
      ;;
    d)
      is_do_not_merge_A="true"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

# parse arguments
# ---------------
# shift arguments away so that "$@" only contains positional arguments
shift "$((OPTIND-1))"

# get input parameters
bed_file_A=${1:-/dev/null}
fasta_fai=${2:-/dev/null}
slop_bp=${3:-1}
bin_size_bp=${4:-1}
direction=${5:-}
op_dir=${6:-}
A_score_filter=${7:-NA}

# perform input checks
# --------------------
# check that both -s and -S are not set
if [ "$num_s_set" -gt 1 ]; then
  >&2 echo "ERROR: only one of -s or -S can be set"
  exit 1
fi

# check that fasta index file exists
if [ ! -f "$fasta_fai" ]; then
  >&2 echo "ERROR: $fasta_fai does not exist"
  exit 1
fi

# check that slop bp is a multiple of bin size bp and that both are positive
if [ $((slop_bp % bin_size_bp)) -ne 0 ]; then
  >&2 echo "ERROR: slop_bp must be a multiple of bin_size_bp"
  exit 1
fi
if [ "$slop_bp" -le 0 ] || [ "$bin_size_bp" -le 0 ]; then
  >&2 echo "ERROR: slop_bp and bin_size_bp must be positive"
  exit 1
fi

# check that A score filter is a number if it is not NA
if [ "$A_score_filter" != "NA" ] && [ ! "$A_score_filter" -eq "$A_score_filter" ] 2>/dev/null; then
  >&2 echo "ERROR: A_score_filter must be a number"
  exit 1
fi

# check that bed file A exists and is a valid bed file
if [ ! -f "$bed_file_A" ]; then
  >&2 echo "ERROR: $bed_file_A does not exist"
  exit 1
fi

if [ ! "$(< "$bed_file_A" python validate_bed_against_fai.py "$fasta_fai")" == "valid" ]; then
  >&2 echo "Error: bed file $bed_file_A must have contigs only in the fasta index file."
  exit 1;
fi

if [ "$intersect_option" == "-s" ] || [ "$intersect_option" == "-S" ] || [ "$A_score_filter" != "NA" ]; then
  if [ ! "$(< "$bed_file_A" python validate_bed_format.py --allow-float-score --six-columns --no-dot-strand)" == "valid" ]; then
    >&2 echo "Error: bed file $bed_file_A is invalid."
    exit 1;
  fi
else
  if [ ! "$(< "$bed_file_A" python validate_bed_format.py --allow-float-score)" == "valid" ]; then
    >&2 echo "Error: bed file $bed_file_A is invalid."
    exit 1;
  fi
fi

# read bed file B from stdin, sort, and send into a temporary file
bed_file_B="$tmpDir"/bed_file_B.bed
cat |  grep -E -v '^#browser|^#track|^#' | sort -k 1,1 -k2,2n > "$bed_file_B"

# check that bed file B is a valid bed file
if [ ! "$(< "$bed_file_B" python validate_bed_format.py --six-columns --no-dot-strand --allow-float-score)" == "valid" ]; then
  >&2 echo "Error: bed file $bed_file_B is invalid."
  exit 1;
fi
if [ ! "$(< "$bed_file_B" python validate_bed_against_fai.py "$fasta_fai")" == "valid" ]; then
  >&2 echo "Error: bed file $bed_file_B must have contigs only in the fasta index file."
  exit 1;
fi

# check that direction is one of the allowed values
if [ "$direction" != "head-to-tail" ] && [ "$direction" != "tail-to-head" ] && \
  [ "$direction" != "increasing-ref" ] && [ "$direction" != "decreasing-ref" ]; then
  >&2 echo "ERROR: direction must be one of 'head-to-tail', 'tail-to-head', 'increasing-ref', 'decreasing-ref'"
  exit 1;
fi

# make the output directory if it does not exist
mkdir -p "$op_dir"

# grow the intervals of bed file B by slop_bp
# -------------------------------------------
< "$bed_file_B" python grow_and_split_bed_file_by_column_value_and_length.py --output-prefix "$op_dir"/bin \
    --fai-file "$fasta_fai" --align-by "$direction" --grow-region-num-bases "$slop_bp" \
    --length-split-num-bases "$bin_size_bp" > "$tmpDir"/split.json

# merge the intervals of bed file A
# ----------------------------------
# we will merge the intervals of bed file A so that overlapping or book-ended intervals are treated as one.
proc_bed_file_A="$tmpDir"/proc_bed_file_A.bed
< "$bed_file_A" sort -k1,1 -k2,2n | awk -v a_score="$A_score_filter" '$3 > $2 && (a_score == "NA" || $5 >= a_score )' |\
  {
    if [ "$is_do_not_merge_A" == "false" ]; then
      if [ "$intersect_option" == "-s" ] || [ "$intersect_option" == "-S" ]; then
        bedtools merge -i stdin -s -c 6 -o distinct |\
          awk -v OFS="\t" '{print $1,$2,$3,"a",1000,$4}'
      else
        bedtools merge -i stdin |\
          awk -v OFS="\t" '{print $1,$2,$3}'
      fi
    else
      cat
    fi
  } | sort -k 1,1 -k2,2n > "$proc_bed_file_A"

# count the number of intervals of bed file A that overlap with each bin and output to a tab-separated file and stdout
# --------------------------------------------------------------------------------------------------------------------
# output header to the tab-separated file
{
  echo "# Number of entries in A file (after filtration/merging depending on options): $(wc -l < "$proc_bed_file_A")"
  echo -e "filename\tdistance\tcount\tbp_in_bin"
} | insert_calling_script_header "$@" > "$op_dir"/counts_per_bin.tsv

# loop over the bin files
for bin_file in $(jq -r '.[].bed_file' "$tmpDir"/split.json ); do

  # if bin file does not contain the string "LS" then it is not the split file, so skip it
  if [[ ! "$(basename "$bin_file")" == *"LS"* ]]; then
    continue
  fi

  # count the number of intervals of bed file A that overlap with the bin
  count=$(awk '$3 > $2' "$bin_file" |\
    bedtools intersect -u $intersect_option -a "$proc_bed_file_A" -b - | wc -l)

  # count the number of base pairs in the bin
  bp_in_bin=$(awk '$3 > $2' "$bin_file" | bedtools merge -i stdin -s -c 6 -o distinct |\
    awk -v sum=0 '{sum+=$3-$2} END{print sum}')

  # get the distance
  if [ "$direction" == "head-to-tail" ] || [ "$direction" == "tail-to-head" ]; then
    distance_option="b"
  elif [ "$direction" == "increasing-ref" ] || [ "$direction" == "decreasing-ref" ]; then
    distance_option="ref"
  fi

  bin_distance=$( grep -v '^#' "$bin_file" | awk -v OFS="\t" '{if($3>$2){b=($2+$3)/2;print $1, b, b, $4, $5, $6}}' |\
    bedtools closest -a - -b "$bed_file_B" -D "$distance_option" | awk '$1 == $7 {print $NF}' |\
      mlr --tsv --implicit-csv-header --headerless-csv-output stats1 -a mean -f 1)

  if [ "$direction" == "increasing-ref" ] || [ "$direction" == "head-to-tail" ]; then
    bin_distance=$(echo "$bin_distance" | awk '{print -$1}')
  fi

  # output the filename, distance, and count to a tab-separated file
  echo -e "$(basename "$bin_file")\t$bin_distance\t$count\t$bp_in_bin" >> "$op_dir"/counts_per_bin.tsv
done

# remove temporary directory
rm -r "$tmpDir"

# output the counts per bin to stdout
cat "$op_dir"/counts_per_bin.tsv