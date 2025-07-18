#!/bin/bash
#SBATCH --mem-per-cpu=10G
#SBATCH -c 10
#SBATCH -p ei-medium
#SBATCH -J bedRegShortPauseReportPerEntry
#SBATCH --mail-type=END,FAIL
#SBATCH --time=23:59:59
#SBATCH --constraint=""

# goal
# -----
# Per entry in a bed file, report pause statistics such as number of pauses, number of expected pauses etc.
# This is a shorter version of the bed_region_pause_report.sh that outputs many other statistics and collects
# them into a pdf and operates over an entire bed file instead of per entry.
# If you want anything further than what is reported here (e.g.: read ids of the pauses corresponding to
# a bed entry of interest), you should use bed_region_pause_report.sh.

# usage
#------
# sbatch bed_region_short_pause_report_per_entry.sh bed_file pause_file fasta_file \
#   pause_sensitivity_bedgraphs_prefix_optional relative_direction_optional restrict_fork_direction_optional
# bed_file: bed file with regions of interest, must be in the six column format.
#           NOTE: this script does not check for overlaps within the bed file etc. and just reports some values
#                 per entry. It is up to the user to ensure that the bed file makes sense.
# pause_file: Pause information in our tab-separated value format. Must contain the columns detectIndex, pauseSite at
#             minimum. For more details on the format, consult the script convert_pause_file_to_bed.py.
# fasta_file: Reference genome in fasta format.
# pause_sensitivity_bedgraphs_prefix_optional: (optional, default "" i.e. unused)
#                                              Prefix of the bedgraphs that contain the pause sensitivity information.
#                                              Our pause pipeline produces a sophisticated null hypothesis for how many
#                                              pauses are expected vs reference position. This is stored in several
#                                              bedgraphs named /A/B/C/D.suffix.bedgraph where suffix can be "all",
#                                             "left", "right", "left.minus" etc. For this example, set this option to
#                                             /A/B/C/D . If this option is not provided, the script will not output
#                                             the expected number of pauses according to the null hypothesis.
# relative_direction_optional: (optional, default "all") type of relative directions to allow,
#                               set to "all", "co-directional" or "head-on". The adjectives refer to the relative
#                               orientation between a fork and a bed entry.
# restrict_fork_direction_optional: (default "" i.e. no restriction) restrict fork direction to "L"/"R"/"lead"/"lag"
#                                    i.e. only use left, right, leading or lagging forks respectively in any
#                                    calculation. We know that in reality forks perform both leading and lagging strand
#                                    synthesis. As our data is _single-molecule_, we _can_ classify forks as leading
#                                    or lagging. Default is no restriction.

# outputs
# -------
# Input bed file is output with five additional columns to standard output:
# As, Ts, pause_count, mean_pause_duration, pause_sensitivity
# NOTE: although our method only calls pauses at Ts, we count pauses on either strand that overlap with a feature
#       as being caused by that feature. So both As and Ts are counted and reported here in separate columns.
#       Exceptions are if you are interested in a combination of a relative direction and a restrict direction
#       like head on-lead. In this case, only the A column or the T column matter. We leave it to the user
#       to think about this downstream of this script.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages
source load_package.sh -samtools -jq -python -miller -bedtools

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

# assign arguments to variables
bed_file="$1"
pause_file="$2"
fasta_file="$3"
pause_sensitivity_bedgraphs_prefix="${4:-}"
relative_direction="${5:-all}"
restrict_fork_direction="${6:-}"

# check that the correct number of arguments were provided
if [ "$#" -lt 3 ]; then
  echo >&2 "ERROR: incorrect number of arguments provided"
  echo >&2 "USAGE: sbatch <script_name.sh> bed_file pause_file fasta_file \ "
  echo >&2 "         pause_sensitivity_bedgraphs_prefix_optional relative_direction_optional \ "
  echo >&2 "         restrict_fork_direction_optional"
  echo >&2 "For more details on what these parameters mean, see the comments in the script."
  echo >&2 "The suffix _optional means that the parameter is optional."
  exit 1
fi

# check that the input files exist
if [ ! -f "$bed_file" ] || [ ! -f "$pause_file" ] || [ ! -f "$fasta_file" ] || [ ! -f "$fasta_file.fai" ]; then
  echo >&2 "ERROR: One of the input files or an associated file (like .fai) does not exist."
  exit 1
fi

# if pause_sensitivity_bedgraphs_prefix is provided, check that the files exist
if [ ! "$pause_sensitivity_bedgraphs_prefix" == "" ]; then
  for suffix in all left right left.minus right.minus left.plus right.plus; do
    if [ ! -f "$pause_sensitivity_bedgraphs_prefix.$suffix.bedgraph" ]; then
      echo >&2 "ERROR: pause_sensitivity_bedgraphs_prefix is provided but one or more files do not exist."
      exit 1
    fi
  done
fi

# check that the input files are valid
if [ ! "$(< "$bed_file" python validate_bed_format.py --six-columns --allow-float-score)" == "valid" ]; then
  >&2 echo "Error: bed file is not valid."
  exit 1;
fi

# bed file must have contigs only in the fasta index file
if [ ! "$(< "$bed_file" python validate_bed_against_fai.py "$fasta_file".fai)" == "valid" ]; then
  >&2 echo "Error: bed file must have contigs only in the fasta index file."
  exit 1;
fi

# check that the pause file is valid
if [ ! "$(< "$pause_file" python validate_pause_format.py)" == "valid" ]; then
  >&2 echo "Error: pause file is not valid."
  exit 1;
fi

# create a temporary bed file with only the true pauses by discarding any pause where any keep column is not True.
# The columns are: contig, start, end, read id, score set to 1, -/+ depending on whether fork is L/R,
# -/+ depending on alignment strand, and pause duration.
tmp_true_pause_bed_file=$(mktemp "$tmpDir"/XXXXXXXXXX)
tmp_true_pause_bed_file_step_1=$(mktemp "$tmpDir"/XXXXXXXXXX)
< "$pause_file" python convert_pause_file_to_bed.py --discardUnkeptPauses --outputDurationAsScore \
    --outputBed6Plus1WherePlus1isAlignStrand --LRtoPlusMinus |\
      grep -E -v '^browser|^track|^#' |\
        awk 'BEGIN{OFS="\t"}{print $1, $2, $3, $4, 1, $6, $7, $5}' |\
          bedtools sort -g "$fasta_file.fai" > "$tmp_true_pause_bed_file_step_1"

# create a temporary bed file with using the sensitivity bedgraphs.
# The format of the bed file is similar to the temporary pause bed file above.
tmp_pause_sensitivity_bed_file=$(mktemp "$tmpDir"/XXXXXXXXXX)
tmp_pause_sensitivity_bed_file_step_1=$(mktemp "$tmpDir"/XXXXXXXXXX)
if [ ! "$pause_sensitivity_bedgraphs_prefix" == "" ]; then
  suffix_list=(left.minus right.minus left.plus right.plus)
  direction_list=(- + - +)
  strand_list=(- - + +)

  # Check that the sensitivity bedgraph file is valid.
  # We want to check that the sensitivity file is "genome covering" i.e. there are windows that cover
  # every base in the genome (we tolerate some missing bases) and that there is no overlap between windows.
  # Our later calculations rely on there being a value for every base in the genome.
  for suffix in "${suffix_list[@]}"; do

    tmp_validity_check=$(mktemp "$tmpDir"/XXXXXXXXXX)
    < "$pause_sensitivity_bedgraphs_prefix.$suffix.bedgraph" grep -E -v '^browser|^track|^#|^chrM' |\
        awk  -v IFS=" " -v OFS="\t" '{print $1, $2, $3}' |\
          python validate_bed_against_fai.py --check-genome-cov-to-tol 100 \
            --check-no-overlap --exclude-contigs chrM "$fasta_file.fai" > "$tmp_validity_check"

    if [ ! "$(cat "$tmp_validity_check")" == "valid" ]; then
      >&2 echo "Error: sensitivity file is not valid. Please check the script to see what we mean."
      exit 1;
    fi
  done

  # We convert the bedgraph file into a bed file - the first four columns are contig, start, end, and a name
  # of "blank". The score column is the bedgraph signal normalized by the window size.
  # The strand is -/+ depending on whether the bedgraph contained left/right fork data - this is useful
  # for us to do head-on/co-directional calculations.
  # The next column is +/- depending on whether the bedgraph corresponded to plus/minus strand of the reference genome.
  {
    for i in {0..3}; do

      suffix="${suffix_list[$i]}"
      direction="${direction_list[$i]}"
      strand="${strand_list[$i]}"

      < "$pause_sensitivity_bedgraphs_prefix.$suffix.bedgraph" grep -E -v '^browser|^track|^#' |\
        awk -v IFS=" " -v OFS="\t" -v d="$direction" -v s="$strand" '{print $1, $2, $3, "blank", $4/($3-$2), d, s}'

    done
  } | bedtools sort -g "$fasta_file.fai" > "$tmp_pause_sensitivity_bed_file_step_1"
else
  : # do nothing
fi

# FLAG: LEAD LAG DISTINCTION
if [ "$restrict_fork_direction" == "L" ]; then
  awk '$6=="-"' "$tmp_true_pause_bed_file_step_1" > "$tmp_true_pause_bed_file"
  awk '$6=="-"' "$tmp_pause_sensitivity_bed_file_step_1" > "$tmp_pause_sensitivity_bed_file"
elif [ "$restrict_fork_direction" == "R" ]; then
  awk '$6=="+"' "$tmp_true_pause_bed_file_step_1" > "$tmp_true_pause_bed_file"
  awk '$6=="+"' "$tmp_pause_sensitivity_bed_file_step_1" > "$tmp_pause_sensitivity_bed_file"
elif [ "$restrict_fork_direction" == "lead" ]; then
  awk '$6==$7 && ($6=="+" || $6=="-")' "$tmp_true_pause_bed_file_step_1" > "$tmp_true_pause_bed_file"
  awk '$6==$7 && ($6=="+" || $6=="-")' "$tmp_pause_sensitivity_bed_file_step_1" > "$tmp_pause_sensitivity_bed_file"
elif [ "$restrict_fork_direction" == "lag" ]; then
  awk '$6!=$7 && ($6=="+" || $6=="-")' "$tmp_true_pause_bed_file_step_1" > "$tmp_true_pause_bed_file"
  awk '$6!=$7 && ($6=="+" || $6=="-")' "$tmp_pause_sensitivity_bed_file_step_1" > "$tmp_pause_sensitivity_bed_file"
elif [ "$restrict_fork_direction" == "" ]; then
   cp "$tmp_true_pause_bed_file_step_1" "$tmp_true_pause_bed_file"
   cp "$tmp_pause_sensitivity_bed_file_step_1" "$tmp_pause_sensitivity_bed_file"
else
  echo >&2 "ERROR: restrict_fork_direction must be one of L, R, lead, lag or empty."
  exit 1
fi

# insert the script header
echo '# last few columns are: n_A, n_T, n_pauses, mean_pause_duration, pause_sensitivity' |\
  insert_calling_script_header "$@"

echo "# NOTE: the A and T columns are output so that the user can calculate a null hypothesis based on "
echo "#       pauses being equally distributed among thymidines."
echo "# NOTE: In general, we associate a pause with a bed file entry if there is same-strand or opposite-strand"
echo "#       overlap, unless the user requests a specific relative direction and a fork combination such as"
echo "#       head-on + lead etc."
echo "#       Because of this, and other reasons like there may be overlap between lines in the bed file, "
echo "#       the user should be careful about adding results across different bed lines. Let us say a pause "
echo "#       is found to overlap with multiple bed lines, then the pause is counted for each bed line in this script."

# find number of self-intersections for informational purposes.
# This is not used elsewhere in the script.
# first, extract first six columns from the bed file
tmp_bed_file=$(mktemp "$tmpDir"/XXXXXXXXXX)
< "$bed_file" grep -E -v '^browser|^track|^#' |\
  sort -k 1,1 -k2,2n |\
  awk 'BEGIN{OFS="\t"}{if($2!=$3){print $1,$2,$3,$4,$5,$6}}' > "$tmp_bed_file"
# shellcheck disable=SC1010
# shellcheck disable=SC2016
twice_number_of_intersections=\
$(bedtools intersect -a "$tmp_bed_file" -b "$tmp_bed_file" -wo |\
   mlr --tsv --skip-comments --implicit-csv-header --headerless-csv-output \
   filter '!($2 == $8 && $3 == $9 && $6 == $12)' then count)
echo "# number_of_intersections_ignoring_strand_and_ignoring_self: $((twice_number_of_intersections/2))"
# shellcheck disable=SC1010
# shellcheck disable=SC2016
twice_number_of_intersections=\
$(bedtools intersect -a "$tmp_bed_file" -b "$tmp_bed_file" -wo -s |\
   mlr --tsv --skip-comments --implicit-csv-header --headerless-csv-output \
   filter '!($2 == $8 && $3 == $9 && $6 == $12)' then count)
echo "# number_of_intersections_retaining_strand_and_ignoring_self: $((twice_number_of_intersections/2))"
echo "# The number of self-intersections is not used elsewhere in the script. It is just for informational purposes."

# count number of columns in the bed file and make sure every line has the same number of columns
column_counts=$(< "$bed_file" grep -E -v '^browser|^track|^#' | awk -F'\t' '{print NF}' | sort -u)
unique_counts=$(echo "$column_counts" | wc -l)
if [ "$unique_counts" -ne 1 ]; then
  echo >&2 "ERROR: bed file has different number of columns in different lines."
  exit 1;
else
  echo "# number_of_columns_in_bed_file: $column_counts"
fi

# set strand option for intersects depending on the relative direction input
if [ "$relative_direction" == "all" ]; then
  strand_option=""
elif [ "$relative_direction" == "co-directional" ]; then
  strand_option="-s"
elif [ "$relative_direction" == "head-on" ]; then
  strand_option="-S"
else
  echo >&2 "ERROR: relative_direction must be one of all, co-directional, or head-on."
  exit 1
fi

# calculate number of lines with entries of zero size and non-zero size respectively and print them
num_non_zero_entries=$(wc -l < "$tmp_bed_file")
num_total_entries=$(grep -c -E -v '^browser|^track|^#' "$bed_file")
num_zero_entries=$((num_total_entries - num_non_zero_entries))
echo "# number of entries with zero size: $num_zero_entries"
echo "# number of entries with non zero size: $num_non_zero_entries"

# perform intersection to calculate all quantities except pause sensitivity
# In order, we do:
# - remove comments, header, and any bed entry with zero size.
# - add a column with the line number, this is helpful to remember the original order of the lines in the bed file
# - sort the bed file
# - count A,G,C,T,N bases in the bed file and retain only A and T counts,
# - gather intersections with the pause file and sum its 5th and 8th columns. If you recollect, column 5 is identically
#   1, and column 8 is the pause duration. So the sum gives total number of pauses and total pause duration.
#   We do not need to set overlap parameters for bedtools map as the pause bed file has intervals of size 1.
# - Print all columns except the last and replace the last column with the mean pause duration (or NA if no pauses)
if [ "$num_non_zero_entries" -gt 0 ]; then
  tmp_output_bed_file=$(mktemp "$tmpDir"/XXXXXXXXXX)

  < "$bed_file" grep -E -v '^browser|^track|^#' |\
    awk -v OFS="\t" '{if($2 != $3){print $0, NR}}' |\
      bedtools sort -g "$fasta_file.fai" |\
        python count_bases_windows_given_fasta_and_bed.py -i stdin -g "$fasta_file" |\
          awk '{for (i=1;i<NF-4;i++) printf "%s\t",$i; printf "%s\t%s\n", $((NF-4)), $((NF-1)) }' |\
            bedtools map $strand_option -a - -b "$tmp_true_pause_bed_file" -c 5,8 -o sum,sum \
              -g "$fasta_file.fai" -null 0 |\
                awk '{for (i=1;i<=NF-1;i++) printf "%s\t",$i; printf "%s\n", ($((NF-1)) > 0) ? $NF/$((NF-1)) : "NA"}' \
                  > "$tmp_output_bed_file"
fi

# add sensitivity values (if requested),
# restore bed entries of zero size after adding NAs for all measurements reported here,
# sort by original line number,
# and delete the column with the original line number.
{
  if [ "$num_non_zero_entries" -gt 0 ]; then
    if [ "$pause_sensitivity_bedgraphs_prefix" == "" ]; then
      < "$tmp_output_bed_file" awk -v OFS="\t" '{print $0, "NA"}'
    else
      # calculate sensitivity.
      # First, we gather all sensitivity windows that overlap with the bed entry,
      #        then we multiply the sensitivity score per base with the amount of overlap,
      #        and finally we sum the values.
      #        We also get rid of columns corresponding to contig, start, end etc. of the sensitivity bed file.
      # The sensitivity score is the expected number of pauses per base per window according to the null hypothesis.
      # (Earlier in this script, we normalized each window's sensitivity to the window size.)
      bedtools intersect $strand_option -a "$tmp_output_bed_file" -b "$tmp_pause_sensitivity_bed_file" -wo |\
          awk -v OFS="\t" '{print $0, $NF * $((NF - 3))}' |\
            cut -f 1-"$((column_counts + 5)),$((column_counts + 5 + 9))" |\
              bedtools groupby -g 1-"$((column_counts + 5))" -c "$((column_counts + 6))" -o sum
    fi
  fi

  # Output bed entries of zero size
  if [ "$num_zero_entries" -gt 0 ]; then
    < "$bed_file" grep -E -v '^browser|^track|^#' |\
      awk -v OFS="\t" '{if($2 == $3){print $0, NR, "NA", "NA", "NA", "NA", "NA"}}';
  fi

} | sort -k"$((column_counts+1)),$((column_counts+1))"n | cut -f "$((column_counts+1))" --complement

# remove temporary directory
rm -rf "$tmpDir"