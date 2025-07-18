#!/bin/bash

# goal
# -----
# count reads in rDNA region that correspond to parent, nascent, fwd, rev

# usage
#------
# bash rDNA_count_parent_nascent_fwd_rev.sh <parent.bed> <nascent.bed> <pauseFile> <fai_file> <contig> <n_repeat_unit>
# no need to use sbatch as this script is not that computationally intensive
# parent.bed: parental bed file produced by rDNA pipeline
# nascent.bed: nascent bed file produced by rDNA pipeline
# pauseFile: pause file produced by rDNA pipeline, in our usual pause format i.e. tab-separated with column headers
#            and must contain the columns detectIndex, pauseSite
# fai_file: fasta index file for the reference genome
# contig: contig name which contains the repeat unit
# n_repeat_unit: number of repeat units in the contig

# outputs
# -------
# tab-separated format with headers and column names is printed to stdout.
# columns are: category, count

# stop execution if any command fails
set -e

# load git repo labels
source load_git_repo_labels.sh

# load packages 
source load_package.sh -python -miller

# load configuration
source config.sh

# make a few temporary files
tmp_dir=$(mktemp -d)
temp_parent_nascent_combine=$(mktemp --tmpdir="$tmp_dir")
temp_pause_bed=$(mktemp --tmpdir="$tmp_dir")
temp_file_parent_readids=$(mktemp --tmpdir="$tmp_dir")
temp_file_nascent_readids=$(mktemp --tmpdir="$tmp_dir")
temp_file_pn_common_readids=$(mktemp --tmpdir="$tmp_dir")
temp_file_pause_readids=$(mktemp --tmpdir="$tmp_dir")

# function to enforce periodicity of the pause sites and restrict to six columns
enforce_periodicity_bed_6(){
  n_repeat_unit_bases=$1
  awk -F'\t' -v n="$n_repeat_unit_bases" -v OFS='\t' '{print $1, $2 - int($2/n)*n, $3 - int($3/n)*n, $4, $5, $6}'
}

# function to set calling script information
insert_calling_script_header() {
  sed '1i'\
'# from commit '"${COMMITSTR:-NA}"' generated at '"${TIMENOW:-NA}"' by '"${config[name]:-NA}"' <'"${config[email]:-NA}"'>\n'\
'# script: '"$0"'\n'\
'# arguments: '"$*"'\n'\
"# slurm job name: ${SLURM_JOB_NAME:-NA}"
}

# check that the correct number of arguments were provided
if [ "$#" -ne 6 ]; then
    >&2 echo "ERROR: incorrect number of arguments provided"
    >&2 echo "usage: bash rDNA_count_parent_nascent_fwd_rev.sh <parent.bed> <nascent.bed> <pauseFile> <fai_file> <contig> <n_repeat_unit>"
    >&2 echo "parent.bed: parental bed file produced by rDNA pipeline"
    >&2 echo "nascent.bed: nascent bed file produced by rDNA pipeline"
    >&2 echo "pauseFile: pause file produced by rDNA pipeline, in our usual pause format i.e. tab-separated with column headers"
    >&2 echo "           and must contain the columns detectIndex, pauseSite"
    >&2 echo "fai_file: fasta index file for the reference genome"
    >&2 echo "contig: contig name which contains the repeat unit"
    >&2 echo "n_repeat_unit: number of repeat units in the contig"
    exit 1;
fi

# assign arguments to variables
parent_bed=${1:-}
nascent_bed=${2:-}
pause_file=${3:-}
fai_file=${4:-}
contig=${5:-InvalidName}
n_repeat_unit=${6:-0}

# check that the fasta index file exists
if [ ! -f "$fai_file" ]; then
  >&2 echo "ERROR: $fai_file does not exist"
  exit 1;
fi

# check that these are valid bed files and that they exist
for bed_file in "$parent_bed" "$nascent_bed"; do

  if [ ! -f "$bed_file" ]; then
    >&2 echo "ERROR: $bed_file does not exist"
    exit 1;
  fi

  if [ ! "$(< "$bed_file" python validate_bed_format.py --allow-float-score --require-uuid --six-columns)" == "valid" ]; then
      >&2 echo "Error: $bed_file is not in the correct format."
      exit 1;
  fi

  if [ ! "$(< "$bed_file" python validate_bed_against_fai.py "$fai_file" )" == "valid"  ]; then
    >&2 echo "Error: $bed_file does not have valid coordinates."
    exit 1;
  fi

done

# check that pause file exists and is valid
if [ ! -f "$pause_file" ] || [ ! "$(< "$pause_file" python validate_pause_format.py)" == "valid" ]; then
  >&2 echo "ERROR: $pause_file does not exist or is not in the correct format"
  exit 1;
fi

# calculate number of bases in the repeat unit
total_bases_fai=$(< "$fai_file" awk -v contig="$contig" -F'\t' '{if ($1 == contig) {print $2}}')
n_repeat_unit_bases=$((total_bases_fai / n_repeat_unit))
# NOTE: we assume above that the total number of bases in the contig is divisible by the number of repeat units.
# otherwise, we will be off by a few bases.

# combine the bed files into a temporary file
{
  < "$parent_bed" python print_valid_data_bed_lines.py | awk -F'\t' -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, "parent"}'
  < "$nascent_bed" python print_valid_data_bed_lines.py | awk -F'\t' -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, "nascent"}'
} > "$temp_parent_nascent_combine"

# count total number of reads
n_total=$(< "$temp_parent_nascent_combine" wc -l)

# convert the pause file into a temporary file in the bed format
< "$pause_file" python convert_pause_file_to_bed.py --discardUnkeptPauses \
   --LRtoPlusMinus --outputBed6Plus1WherePlus1isAlignStrand |\
     python print_valid_data_bed_lines.py | enforce_periodicity_bed_6 $n_repeat_unit_bases > "$temp_pause_bed"

# count number of pauses
n_pause=$(< "$temp_pause_bed" wc -l)

# calculate number of pauses per bp of repeat unit
n_pause_per_bp_of_repeat=$(echo "scale=6; $n_pause / $n_repeat_unit_bases" | bc)

# count number of pauses on left and right forks and the fraction of right forks pauses
n_pause_left=$(< "$temp_pause_bed" awk -F'\t' '{if ($6 == "-") {print $0}}' | wc -l)
n_pause_right=$(< "$temp_pause_bed" awk -F'\t' '{if ($6 == "+") {print $0}}' | wc -l)
n_pause_right_fraction=$(echo "scale=2; $n_pause_right / ($n_pause_right + $n_pause_left)" | bc)
n_pause_left_per_bp_of_repeat=$(echo "scale=6; $n_pause_left / $n_repeat_unit_bases" | bc)
n_pause_right_per_bp_of_repeat=$(echo "scale=6; $n_pause_right / $n_repeat_unit_bases" | bc)

# calculate statistics for pauses with x coordinate above 8200
# (this is one of the x limits used in our figures where we zoom in near the rRFB).
n_pause_above_8200=$(< "$temp_pause_bed" awk -F'\t' '{if ($2 > 8200) {print $0}}' | wc -l)
n_pause_left_above_8200=$(< "$temp_pause_bed" awk -F'\t' '{if ($6 == "-" && $2 > 8200) {print $0}}' | wc -l)
n_pause_right_above_8200=$(< "$temp_pause_bed" awk -F'\t' '{if ($6 == "+" && $2 > 8200) {print $0}}' | wc -l)

# count pauses between 8755 and 8955 (i.e. 100 bp +- the first Fob1 binding site)
n_pause_between_8755_8955=$(< "$temp_pause_bed" awk -F'\t' '{if ($2 > 8755 && $2 < 8955) {print $0}}' | wc -l)
n_pause_per_bp_between_8755_8955=$(echo "scale=6; $n_pause_between_8755_8955 / 200" | bc)
enrichment_factor_between_8755_8955=$(echo "scale=6; $n_pause_per_bp_between_8755_8955 / $n_pause_per_bp_of_repeat" |bc)

# count pauses from -110 to +2000 on right forks (-110 is the TTS and we are going about 2 kb into the gene)
# we use these statistics in our Fob1D analysis
n_pause_right_between_m110_2000=$(< "$temp_pause_bed" awk -F'\t' -v n="$n_repeat_unit_bases" \
  '{if ($6 == "+" && ($2 > n-110 || $2 < 2000)) {print $0}}' | wc -l)
n_pause_right_per_bp_between_m110_2000=$(echo "scale=6; $n_pause_right_between_m110_2000 / 2110" | bc)
enrichment_factor_right_between_m110_2000=$(echo "scale=6; $n_pause_right_per_bp_between_m110_2000 \
                                              / $n_pause_right_per_bp_of_repeat" | bc)

# count pauses from 4637 to 6747 on left forks (6747 is the TSS and we are going about 2110 bp into the gene,
#  the same length as the -110 to +2000 above)
# we use these statistics in our Fob1D analysis
n_pause_left_between_4637_6747=$(< "$temp_pause_bed" awk -F'\t' -v n="$n_repeat_unit_bases" \
  '{if ($6 == "-" && ($2 > 4637 && $2 < 6747)) {print $0}}' | wc -l)
n_pause_left_per_bp_between_4637_6747=$(echo "scale=6; $n_pause_left_between_4637_6747 / 2000" | bc)
enrichment_factor_left_between_4637_6747=$(echo "scale=6; $n_pause_left_per_bp_between_4637_6747 \
                                              / $n_pause_left_per_bp_of_repeat" | bc)

# count numbers of leading and lagging strand pauses
< "$pause_file" python convert_pause_file_to_bed.py --discardUnkeptPauses \
  --LeadLagToPlusMinus --outputBed6Plus1WherePlus1isAlignStrand |\
    python print_valid_data_bed_lines.py | enforce_periodicity_bed_6 $n_repeat_unit_bases > "$temp_pause_bed"_lead_lag
n_pause_lead=$(< "$temp_pause_bed"_lead_lag awk -F'\t' '{if ($6 == "+") {print $0}}' | wc -l)
n_pause_lag=$(< "$temp_pause_bed"_lead_lag awk -F'\t' '{if ($6 == "-") {print $0}}' | wc -l)
n_pause_lead_above_8200=$(< "$temp_pause_bed"_lead_lag awk -F'\t' '{if ($6 == "+" && $2 > 8200) {print $0}}' | wc -l)
n_pause_lag_above_8200=$(< "$temp_pause_bed"_lead_lag awk -F'\t' '{if ($6 == "-" && $2 > 8200) {print $0}}' | wc -l)

# calculate mean and standard deviation of pause sites above 8200 behind 8855.
# 8855 is the location of the first Fob1 binding site.
# NOTE: behind does not mean we only include pauses behind 8855, but rather we calculate the distance behind 8855.
mean_pause_lead_above_8200_behind_8855=$(< "$temp_pause_bed"_lead_lag \
  awk -F'\t' '{if ($6 == "+" && $2 > 8200) {print 8855-$2}}' |\
    mlr --tsv --implicit-csv-header --headerless-csv-output stats1 -a mean -f 1)
mean_pause_lag_above_8200_behind_8855=$(< "$temp_pause_bed"_lead_lag \
  awk -F'\t' '{if ($6 == "-" && $2 > 8200) {print 8855-$2}}' |\
    mlr --tsv --implicit-csv-header --headerless-csv-output stats1 -a mean -f 1)
sd_pause_lead_above_8200_behind_8855=$(< "$temp_pause_bed"_lead_lag \
  awk -F'\t' '{if ($6 == "+" && $2 > 8200) {print 8855-$2}}' |\
    mlr --tsv --implicit-csv-header --headerless-csv-output stats1 -a stddev -f 1)
sd_pause_lag_above_8200_behind_8855=$(< "$temp_pause_bed"_lead_lag \
  awk -F'\t' '{if ($6 == "-" && $2 > 8200) {print 8855-$2}}' |\
    mlr --tsv --implicit-csv-header --headerless-csv-output stats1 -a stddev -f 1)

# calculate median of all lead, lag pause sites
median_pause_lead=$(< "$temp_pause_bed"_lead_lag \
  awk -F'\t' '{if ($6 == "+") {print 8855-$2}}' |\
    mlr --tsv --implicit-csv-header --headerless-csv-output stats1 -a median -f 1)
median_pause_lag=$(< "$temp_pause_bed"_lead_lag \
  awk -F'\t' '{if ($6 == "-") {print 8855-$2}}' |\
    mlr --tsv --implicit-csv-header --headerless-csv-output stats1 -a median -f 1)

# how many of the reads in the pause file are unique
cut -f4 "$temp_pause_bed" | sort | uniq > "$temp_file_pause_readids"
n_unique_in_pause=$(< "$temp_file_pause_readids" wc -l)

# find if there are any reads that are both parent and nascent
< "$temp_parent_nascent_combine" awk -F'\t' '{if ($7 == "parent") {print $4}}' | sort | uniq > "$temp_file_parent_readids"
< "$temp_parent_nascent_combine" awk -F'\t' '{if ($7 == "nascent") {print $4}}' | sort | uniq > "$temp_file_nascent_readids"
comm -12 "$temp_file_parent_readids" "$temp_file_nascent_readids" > "$temp_file_pn_common_readids"
n_both=$(< "$temp_file_pn_common_readids" wc -l)

# of these reads marked as both parent and nascent, find how many are in the pause file
n_both_in_pause=$(comm -12 "$temp_file_pn_common_readids" "$temp_file_pause_readids" | wc -l)

# count number of unique reads i.e. number of reads where the fourth column is unique
n_unique_parent=$(< "$temp_file_parent_readids" wc -l)
n_unique_nascent=$(< "$temp_file_nascent_readids" wc -l)
n_unique=$((n_unique_parent + n_unique_nascent))

# calculate various categories of lines
n_parent=$(< "$temp_parent_nascent_combine" awk -F'\t' '{if ($7 == "parent") {print $0}}' | wc -l)
n_nascent=$(< "$temp_parent_nascent_combine" awk -F'\t' '{if ($7 == "nascent") {print $0}}' | wc -l)
n_fwd=$(< "$temp_parent_nascent_combine" awk -F'\t' '{if ($6 == "+") {print $0}}' | wc -l)
n_rev=$(< "$temp_parent_nascent_combine" awk -F'\t' '{if ($6 == "-") {print $0}}' | wc -l)

# calculate parent, nascent, fwd, rev, but counting repeat units now
n_parent_repeat=$(< "$temp_parent_nascent_combine" awk -F'\t' -v n="$n_repeat_unit_bases" -v sum=0 \
  '{if ($7 == "parent") {sum += $3 - $2}} END {print int(sum/n)}')
n_nascent_repeat=$(< "$temp_parent_nascent_combine" awk -F'\t' -v n="$n_repeat_unit_bases" -v sum=0 \
  '{if ($7 == "nascent") {sum += $3 - $2}} END {print int(sum/n)}')
n_fwd_repeat=$(< "$temp_parent_nascent_combine" awk -F'\t' -v n="$n_repeat_unit_bases" -v sum=0 \
  '{if ($6 == "+") {sum += $3 - $2}} END {print int(sum/n)}')
n_rev_repeat=$(< "$temp_parent_nascent_combine" awk -F'\t' -v n="$n_repeat_unit_bases" -v sum=0 \
  '{if ($6 == "-") {sum += $3 - $2}} END {print int(sum/n)}')

# calculate total number of repeat units
n_total_repeat=$((n_parent_repeat + n_nascent_repeat))

# calculate pause count normalized to number of repeat units
n_pause_repeat=$(echo "scale=6; $n_pause / $n_total_repeat" | bc)

# remove temporary files
rm -rf "$tmp_dir"

# output measurements
echo -e "# NOTE: behind 8855 means 'distance behind 8855' and not 'filter to include only pauses behind 8855'"
{
  echo -e "repeat_unit_size\t$n_repeat_unit_bases"
  echo -e "parent_count\t$n_parent"
  echo -e "nascent_count\t$n_nascent"
  echo -e "fwd_count\t$n_fwd"
  echo -e "rev_count\t$n_rev"
  echo -e "total_count\t$n_total"
  echo -e "unique_count\t$n_unique"
  echo -e "unique_parent_count\t$n_unique_parent"
  echo -e "unique_nascent_count\t$n_unique_nascent"
  echo -e "unique_reads_marked_as_both_parent_nascent_count\t$n_both"
  echo -e "pause_count\t$n_pause"
  echo -e "pause_count_per_bp_of_repeat\t$n_pause_per_bp_of_repeat"
  echo -e "pause_count_left\t$n_pause_left"
  echo -e "pause_count_right\t$n_pause_right"
  echo -e "pause_count_lead\t$n_pause_lead"
  echo -e "pause_count_lag\t$n_pause_lag"
  echo -e "pause_fraction_right\t$n_pause_right_fraction"
  echo -e "pause_count_above_8200\t$n_pause_above_8200"
  echo -e "pause_count_left_above_8200\t$n_pause_left_above_8200"
  echo -e "pause_count_right_above_8200\t$n_pause_right_above_8200"
  echo -e "pause_count_lead_above_8200\t$n_pause_lead_above_8200"
  echo -e "pause_count_lag_above_8200\t$n_pause_lag_above_8200"
  echo -e "unique_reads_in_pause_count\t$n_unique_in_pause"
  echo -e "unique_reads_marked_as_both_parent_nascent_in_pause_count\t$n_both_in_pause"
  echo -e "parent_count_repeat_unit\t$n_parent_repeat"
  echo -e "nascent_count_repeat_unit\t$n_nascent_repeat"
  echo -e "fwd_count_repeat_unit\t$n_fwd_repeat"
  echo -e "rev_count_repeat_unit\t$n_rev_repeat"
  echo -e "total_count_repeat_unit\t$n_total_repeat"
  echo -e "pause_count_repeat_unit\t$n_pause_repeat"
  echo -e "mean_pause_lead_above_8200_behind_8855\t$mean_pause_lead_above_8200_behind_8855"
  echo -e "mean_pause_lag_above_8200_behind_8855\t$mean_pause_lag_above_8200_behind_8855"
  echo -e "sd_pause_lead_above_8200_behind_8855\t$sd_pause_lead_above_8200_behind_8855"
  echo -e "sd_pause_lag_above_8200_behind_8855\t$sd_pause_lag_above_8200_behind_8855"
  echo -e "pause_count_between_8755_8955\t$n_pause_between_8755_8955"
  echo -e "pause_count_per_bp_between_8755_8955\t$n_pause_per_bp_between_8755_8955"
  echo -e "enrichment_factor_between_8755_8955\t$enrichment_factor_between_8755_8955"
  echo -e "pause_count_right_between_m110_2000\t$n_pause_right_between_m110_2000"
  echo -e "pause_count_right_per_bp_between_m110_2000\t$n_pause_right_per_bp_between_m110_2000"
  echo -e "enrichment_factor_right_between_m110_2000\t$enrichment_factor_right_between_m110_2000"
  echo -e "pause_count_left_between_4637_6747\t$n_pause_left_between_4637_6747"
  echo -e "pause_count_left_per_bp_between_4637_6747\t$n_pause_left_per_bp_between_4637_6747"
  echo -e "enrichment_factor_left_between_4637_6747\t$enrichment_factor_left_between_4637_6747"
  echo -e "median_pause_lead_behind_8855\t$median_pause_lead"
  echo -e "median_pause_lag_behind_8855\t$median_pause_lag"
} | sed '1icategory\tcount' | insert_calling_script_header "$@"