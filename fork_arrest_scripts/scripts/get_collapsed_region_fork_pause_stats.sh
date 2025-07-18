#!/bin/bash

# goal
# -----
# Given regions of interest, forks, and pauses, collapse one feature (region or fork) and plot
# coverage of the other feature(s) (pause or region) on it.

# usage
#------
# bash get_collapsed_region_fork_pause_stats.sh all_fork_file pause_file region_bed output_dir
# can use sbatch if need be but not necessary unless the files are large.
# all_fork_file: a file in our pause format that contains pause information from all forks passing through the ROI.
#                The word 'all' means we note down a fork irrespective of whether it has pauses.
#                Such a file is produced if you run the script bed_region_pause_report.sh for example.
#                This file must be in our standard pause format, which is tab-separated with header and column names,
#                and with columns detectIndex, pauseSite, keep* etc. We are not gonna get into the format details here.
#                See comments on other scripts or the validate_pause_format.py file to learn more.
# pause_file: a file in our pause format that contains all pauses in the ROI.
#             Such a file is produced if you run the script bed_region_pause_report.sh for example.
#             Conceptually, this file is a subset of the all_fork_file with valid pauses in the region.
# region_bed: a bed file containing the regions of interest. File must have at least six columns.
# output_dir: the directory where the output will be written. This directory will be created if it doesn't exist.
#             preferable to use an empty directory.
# NOTE: if you do not want to use one of the files, you can use /dev/null instead of the file name.

# outputs
# -------
# several plots are sent to the output directory.

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -R -python -bedtools -miller

# load configuration
source config.sh

# set temporary directory
tmpDir="${config[scratchDir]:-}"/tmp/temp$(openssl rand -hex 6)
mkdir -p "$tmpDir"

# check that the correct number of arguments were provided
if [ "$#" -lt 4 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash get_collapsed_region_fork_pause_stats.sh all_fork_file pause_file region_bed output_dir"
    >&2 echo "For more details on what these parameters mean, see the comments in the script."
    exit 1;
fi

# assign arguments to variables
all_fork_file=${1:-}
pause_file=${2:-}
region_bed=${3:-}
output_dir=${4:-}

# make output directory if it doesn't exist
mkdir -p "$output_dir"

# check that the input files exist and are valid
if [ ! -f "$all_fork_file" ]; then
  >&2 echo "WARNING: all_fork_file $all_fork_file does not exist, so disabling related plots"
  all_fork_file=/dev/null
elif [ ! "$(< "$all_fork_file" python validate_pause_format.py)" == "valid" ]; then
  >&2 echo "WARNING: pause file $all_fork_file is not valid, so disabling related plots"
  all_fork_file=/dev/null
fi

if [ ! -f "$pause_file" ]; then
  >&2 echo "WARNING: pause_file $pause_file does not exist, so disabling related plots"
  pause_file=/dev/null
elif [ ! "$(< "$pause_file" python validate_pause_format.py)" == "valid" ]; then
  >&2 echo "WARNING: pause_file $pause_file is not valid, so disabling related plots"
  pause_file=/dev/null
fi

if [ ! -f "$region_bed" ]; then
  >&2 echo "ERROR: region_bed $region_bed does not exist, so disabling related plots"
  region_bed=/dev/null
elif [ ! "$(< "$region_bed" python validate_bed_format.py --six-columns --allow-float-score)" == "valid" ]; then
  >&2 echo "Error: bed file $region_bed is not valid, so disabling related plots"
  region_bed=/dev/null
fi

if [ -f "$region_bed" ]; then

  # merge bed file by strand
  # ========================
  tmp_merged_bed_file=$(mktemp -p "$tmpDir" XXXXXX.bed);

  sort -k 1,1 -k2,2n "$region_bed" |\
        bedtools merge -s -i stdin -c 6 -o distinct |\
          awk 'BEGIN{OFS="\t"}{if($2!=$3){print $1, $2, $3, "blank", 1000, $4}}' > "$tmp_merged_bed_file"

  # get max region size
  # ===================
  # loop through the bed file and find the length of the largest interval
  max_region_size=$(awk -v a=0 '{if($3-$2>a) a=$3-$2}END{print a}' "$tmp_merged_bed_file")

fi

if [ -f "$all_fork_file" ]; then

  # get max fork size
  # =================
  # loop through the fork bed file and find the length of the largest interval
  max_fork_size=$(< "$all_fork_file" python convert_pause_file_to_bed.py --outputForks |\
    grep -E -v '^browser|^track|^#' |\
    awk -v a=0 '{if($3-$2>a) a=$3-$2}END{print a}')

  # plot fork length distribution
  # =============================
  # shellcheck disable=SC2016
  # shellcheck disable=SC1010
  < "$all_fork_file" python convert_pause_file_to_bed.py --LRtoPlusMinus --outputForks |\
   mlr --tsv --implicit-csv-header --skip-comments put '$fork_len=$3-$2'\
      then histogram -f fork_len --lo 0 --hi "$max_fork_size" --nbins 100 |\
      Rscript plotting_and_short_analyses/plot_histogram.R fork_len_count\
        "$output_dir"/all_fork_lengths.png  "Fork length (b)" Count

fi

if [ -f "$region_bed" ]; then
  # plot region length distribution
  # =============================
  # shellcheck disable=SC2016
  # shellcheck disable=SC1010
  < "$tmp_merged_bed_file" \
   mlr --tsv --implicit-csv-header --skip-comments put '$region_len=$3-$2'\
      then histogram -f region_len --lo 0 --hi "$max_region_size" --nbins 100 |\
      Rscript plotting_and_short_analyses/plot_histogram.R region_len_count\
        "$output_dir"/all_region_lengths.png  "Region length (b)" Count

   # plot distribution of distance between neighbouring entries on the same strand
    # ==============================================================================
    # shellcheck disable=SC1010
    bedtools closest -io -d -s -t first -a "$tmp_merged_bed_file" -b "$tmp_merged_bed_file" |\
      awk -F"\t" '$1 == $7 {print $NF}'  |\
       mlr --tsv --implicit-csv-header histogram -f 1 --lo 0 --hi 10000 --nbins 200 \
        then rename 1_count,count |\
        Rscript plotting_and_short_analyses/plot_histogram.R count\
          "$output_dir"/all_region_neighbouring_distances.png  "Region length (b)" Count
fi

if [ -f "$region_bed" ] && [ -f "$pause_file" ]; then

  # plot pause profile on subsets by relative direction (all, head-on, codirectional) and fork direction
  # ====================================================================================================
  # (all, lead, lag) (9 plots)
  # ==========================

  tmp_collapse_pause_LR_to_PlusMinus=$(mktemp -p "$tmpDir" XXXXXX);
  tmp_collapse_lead_pause_LR_to_PlusMinus=$(mktemp -p "$tmpDir" XXXXXX);
  tmp_collapse_lag_pause_LR_to_PlusMinus=$(mktemp -p "$tmpDir" XXXXXX);

  < "$pause_file" python convert_pause_file_to_bed.py --LRtoPlusMinus |\
   python calculate_collapsed_bed_coverage_by_bed_file.py "$tmp_merged_bed_file" stdin \
     > "$tmp_collapse_pause_LR_to_PlusMinus"

  # FLAG: LEAD LAG DISTINCTION
  # shellcheck disable=SC2016
  mlr --tsv --skip-comments filter '$detectIndex =~ "_fwd_R_" || $detectIndex =~ "_rev_L_"' "$pause_file" |\
   python convert_pause_file_to_bed.py --LRtoPlusMinus |\
   python calculate_collapsed_bed_coverage_by_bed_file.py "$tmp_merged_bed_file" stdin \
     > "$tmp_collapse_lead_pause_LR_to_PlusMinus"

  # FLAG: LEAD LAG DISTINCTION
  # shellcheck disable=SC2016
  mlr --tsv --skip-comments filter '$detectIndex =~ "_fwd_L_" || $detectIndex =~ "_rev_R_"' "$pause_file" |\
   python convert_pause_file_to_bed.py --LRtoPlusMinus |\
   python calculate_collapsed_bed_coverage_by_bed_file.py "$tmp_merged_bed_file" stdin \
     > "$tmp_collapse_lag_pause_LR_to_PlusMinus"

  suffixDirn=("" "_lead" "_lag")
  input_file_list=("$tmp_collapse_pause_LR_to_PlusMinus" \
    "$tmp_collapse_lead_pause_LR_to_PlusMinus" \
    "$tmp_collapse_lag_pause_LR_to_PlusMinus" \
  )

  suffixRelDirn=("" "_codirectional" "_headon")
  filtration_list=("no_intersection" "same_strand" "opposite_strand")
  inversion_list=("-x" "" "" )

  for dirnCount in {0..2}; do
    for relDirnCount in {0..2}; do
      output_file="$output_dir"/pauses_on_collapsed_region"${suffixDirn[$dirnCount]}""${suffixRelDirn[$relDirnCount]}".png
      output_tsv="$output_dir"/pauses_on_collapsed_region"${suffixDirn[$dirnCount]}""${suffixRelDirn[$relDirnCount]}".tsv
      # shellcheck disable=SC2016
      # shellcheck disable=SC1010
      # shellcheck disable=SC2086
      < "${input_file_list[$dirnCount]}" \
       mlr --tsv --implicit-csv-header --skip-comments put '$abs_start=abs($1)' \
          then rename 3,type then cut -o -f abs_start,type then filter ${inversion_list[$relDirnCount]} '$type=="'"${filtration_list[$relDirnCount]}"'"' |\
          tee "$output_tsv" |\
          mlr --tsv histogram -f abs_start --lo 0 --hi "$max_region_size" --nbins 10 |\
          Rscript plotting_and_short_analyses/plot_histogram.R abs_start_count \
            "${output_file}" "Collapsed Pause profile (b)" \
            Count 0,"$max_region_size"
     done
  done
fi

if [ -f "$region_bed" ] && [ -f "$all_fork_file" ]; then

  # perform calculations to plot region profile on collapsed forks
  # ==============================================================
  tmp_region_on_collapsed_forks=$(mktemp -p "$tmpDir" XXXXXX.tsv);

  < "$all_fork_file" python convert_pause_file_to_bed.py --LRtoPlusMinus --outputForks |\
   python calculate_collapsed_bed_coverage_by_bed_file.py stdin "$tmp_merged_bed_file" > \
   "$tmp_region_on_collapsed_forks"

  # these files are used to store calculations
  tmp_region_coverage_on_collapsed_forks_cd=$(mktemp -p "$tmpDir" XXXXXX.tsv);
  tmp_region_coverage_on_collapsed_forks_ho=$(mktemp -p "$tmpDir" XXXXXX.tsv);
  tmp_fork_coverage_on_collapsed_forks_cd=$(mktemp -p "$tmpDir" XXXXXX.tsv);
  tmp_fork_coverage_on_collapsed_forks_ho=$(mktemp -p "$tmpDir" XXXXXX.tsv);

  # plot region profile on collapsed forks (head-on)
  # =================================================
  tmp_fai_file=$(mktemp -p "$tmpDir" XXXXXX.fai);
  echo -e "chr_dummy\t$max_fork_size" > "$tmp_fai_file";

  # shellcheck disable=SC2016
  # shellcheck disable=SC1010
  < "$tmp_region_on_collapsed_forks"\
     mlr --tsv --implicit-csv-header --headerless-csv-output --skip-comments filter '$3=="opposite_strand"'\
          then put '$contig = "chr_dummy"; $start = min(abs($1),abs($2)); $end = max(abs($1),abs($2));'\
          then cut -o -f contig,start,end |\
          sort -k 1,1 -k2,2n  |\
          bedtools genomecov -i stdin -g "$tmp_fai_file" -d |\
          tee "$tmp_region_coverage_on_collapsed_forks_ho" |\
          sed '1icontig\tstart\tcoverage' |\
          Rscript plotting_and_short_analyses/plot_scatter.R start coverage\
            "$output_dir"/region_on_collapsed_forks_headon.png "Region profile on collapsed forks (b)"\
            "Coverage" 0,auto 0,auto 0.2 1 0 1000


  # plot region profile on collapsed forks (co-directional)
  # =======================================================

  # shellcheck disable=SC2016
  # shellcheck disable=SC1010
  < "$tmp_region_on_collapsed_forks"\
     mlr --tsv --implicit-csv-header --headerless-csv-output --skip-comments filter '$3=="same_strand"'\
          then put '$contig = "chr_dummy"; $start = min(abs($1),abs($2)); $end = max(abs($1),abs($2));'\
          then cut -o -f contig,start,end |\
          sort -k 1,1 -k2,2n  |\
          bedtools genomecov -i stdin -g "$tmp_fai_file" -d |\
          tee "$tmp_region_coverage_on_collapsed_forks_cd" |\
          sed '1icontig\tstart\tcoverage' |\
          Rscript plotting_and_short_analyses/plot_scatter.R start coverage\
            "$output_dir"/region_on_collapsed_forks_codirectional.png "Region profile on collapsed forks (b)"\
            "Coverage" 0,auto 0,auto 0.2 1 0 1000

  # plot fork profile on collapsed forks (head-on)
  # =================================================

  # shellcheck disable=SC2016
  # shellcheck disable=SC1010
  < "$tmp_region_on_collapsed_forks"\
     mlr --tsv --implicit-csv-header --headerless-csv-output --skip-comments filter '$3=="opposite_strand"'\
          then put '$contig = "chr_dummy"; $start = 0; $end = $4;'\
          then cut -o -f contig,start,end |\
          sort -k 1,1 -k2,2n  |\
          bedtools genomecov -i stdin -g "$tmp_fai_file" -d |\
          tee "$tmp_fork_coverage_on_collapsed_forks_ho" |\
          sed '1icontig\tstart\tcoverage' |\
          Rscript plotting_and_short_analyses/plot_scatter.R start coverage\
            "$output_dir"/fork_on_collapsed_forks_headon.png "Fork profile on collapsed forks (b)"\
            "Coverage" 0,auto 0,auto 1


  # plot fork profile on collapsed forks (co-directional)
  # =======================================================

  # shellcheck disable=SC2016
  # shellcheck disable=SC1010
  < "$tmp_region_on_collapsed_forks"\
     mlr --tsv --implicit-csv-header --headerless-csv-output --skip-comments filter '$3=="same_strand"'\
          then put '$contig = "chr_dummy"; $start = 0; $end = $4;'\
          then cut -o -f contig,start,end |\
          sort -k 1,1 -k2,2n  |\
          bedtools genomecov -i stdin -g "$tmp_fai_file" -d |\
          tee "$tmp_fork_coverage_on_collapsed_forks_cd" |\
          sed '1icontig\tstart\tcoverage' |\
          Rscript plotting_and_short_analyses/plot_scatter.R start coverage\
            "$output_dir"/fork_on_collapsed_forks_codirectional.png "Fork profile on collapsed forks (b)"\
            "Coverage" 0,auto 0,auto 1

  # plot normalized region profile on collapsed forks (head-on) (normalized to fork count)
  # ======================================================================================

  tmp_region_coverage_bedgraph_ho=$(mktemp -p "$tmpDir" XXXXXX.bedgraph);
  tmp_fork_coverage_bedgraph_ho=$(mktemp -p "$tmpDir" XXXXXX.bedgraph);

  < "$tmp_region_coverage_on_collapsed_forks_ho" awk 'BEGIN{OFS=" "}{print $1, $2 - 1, $2, $3}' >\
    "$tmp_region_coverage_bedgraph_ho"
  < "$tmp_fork_coverage_on_collapsed_forks_ho" awk 'BEGIN{OFS=" "}{print $1, $2 - 1, $2, $3}' >\
    "$tmp_fork_coverage_bedgraph_ho"

  python divide_bedgraph_values.py "$tmp_region_coverage_bedgraph_ho" "$tmp_fork_coverage_bedgraph_ho"\
    remove_zero_by_zero |\
    sed '1icontig start end coverage' |\
      Rscript plotting_and_short_analyses/plot_scatter.R start coverage\
        "$output_dir"/region_normalized_on_collapsed_forks_headon.png "Region profile on collapsed forks (b)"\
        "Normalized coverage" 0,auto 0,auto 1 1 0 1000


  # plot normalized region profile on collapsed forks (co-directional) (normalized to fork count)
  # =============================================================================================

  tmp_region_coverage_bedgraph_cd=$(mktemp -p "$tmpDir" XXXXXX.bedgraph);
  tmp_fork_coverage_bedgraph_cd=$(mktemp -p "$tmpDir" XXXXXX.bedgraph);

  < "$tmp_region_coverage_on_collapsed_forks_cd" awk 'BEGIN{OFS=" "}{print $1, $2 - 1, $2, $3}' >\
    "$tmp_region_coverage_bedgraph_cd"
  < "$tmp_fork_coverage_on_collapsed_forks_cd" awk 'BEGIN{OFS=" "}{print $1, $2 - 1, $2, $3}' >\
    "$tmp_fork_coverage_bedgraph_cd"

  python divide_bedgraph_values.py "$tmp_region_coverage_bedgraph_cd" "$tmp_fork_coverage_bedgraph_cd"\
    remove_zero_by_zero |\
    sed '1icontig start end coverage' |\
      Rscript plotting_and_short_analyses/plot_scatter.R start coverage\
        "$output_dir"/region_normalized_on_collapsed_forks_codirectional.png "Region profile on collapsed forks (b)"\
        "Normalized coverage" 0,auto 0,auto 1 1 0 1000

  # plot fork profile on collapsed region (head-on)
  # =================================================

  tmp_fai_file=$(mktemp -p "$tmpDir" XXXXXX.fai);
  echo -e "chr_dummy\t$((max_region_size + 1))" > "$tmp_fai_file";
  # shellcheck disable=SC2016
  # shellcheck disable=SC1010
  < "$all_fork_file" python convert_pause_file_to_bed.py --LRtoPlusMinus --outputForks |\
   python calculate_collapsed_bed_coverage_by_bed_file.py "$tmp_merged_bed_file" stdin |\
     mlr --tsv --implicit-csv-header --headerless-csv-output --skip-comments filter '$3=="opposite_strand"'\
          then put '$contig = "chr_dummy"; $start = min(abs($1),abs($2)); $end = max(abs($1),abs($2));'\
          then cut -o -f contig,start,end |\
          sort -k 1,1 -k2,2n  |\
          bedtools genomecov -i stdin -g "$tmp_fai_file" -d |\
          sed '1icontig\tstart\tcoverage' |\
          Rscript plotting_and_short_analyses/plot_scatter.R start coverage\
            "$output_dir"/fork_on_collapsed_region_headon.png "Fork profile on collapsed region (b)"\
            "Coverage" 0,auto 0,auto 1


  # plot fork profile on collapsed region (co-directional)
  # =======================================================

  tmp_fai_file=$(mktemp -p "$tmpDir" XXXXXX.fai);
  echo -e "chr_dummy\t$((max_region_size + 1))" > "$tmp_fai_file";
  # shellcheck disable=SC2016
  # shellcheck disable=SC1010
  < "$all_fork_file" python convert_pause_file_to_bed.py --LRtoPlusMinus --outputForks |\
   python calculate_collapsed_bed_coverage_by_bed_file.py "$tmp_merged_bed_file" stdin |\
     mlr --tsv --implicit-csv-header --headerless-csv-output --skip-comments filter '$3=="same_strand"'\
          then put '$contig = "chr_dummy"; $start = min(abs($1),abs($2)); $end = max(abs($1),abs($2));'\
          then cut -o -f contig,start,end |\
          sort -k 1,1 -k2,2n  |\
          bedtools genomecov -i stdin -g "$tmp_fai_file" -d |\
          sed '1icontig\tstart\tcoverage' |\
          Rscript plotting_and_short_analyses/plot_scatter.R start coverage\
            "$output_dir"/fork_on_collapsed_region_codirectional.png "Fork profile on collapsed region (b)"\
            "Coverage" 0,auto 0,auto 1

fi

# remove temporary directory
rm -rf "$tmpDir"