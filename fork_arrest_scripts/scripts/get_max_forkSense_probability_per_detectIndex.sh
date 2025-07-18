#!/usr/bin/env bash
#SBATCH --mem=200G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J maxForkProbPerDI
#SBATCH --mail-type=END,FAIL
#SBATCH --constraint=""
#SBATCH --time=23:59:59

# goal of program
# ================
# Print maximum fork sense probability associated with each detectIndex

# program usage
# =============
# sbatch <program_name>.sh detectIndex_file mod_bam_left mod_bam_right
# can use bash in place of sbatch
# detect_index_file: tab-separated file with headers, one required column detectIndex, other columns ignored.
#                     detectIndex must have the format 'readID_contig_start_end_orientation_dirnFork_startFork_endFork'
#                     orientation=fwd/rev, dirnFork=L/R.
#                     start, end, startFork, endFork are coordinates on the reference genome.
#                     startFork always lesser than endFork irrespective of fork direction.
# mod_bam_left: bam file with left fork probability at each thymidine.
# mod_bam_right: bam file with right fork probability at each thymidine.

# load variables from the command line
pause_file=$1;
mod_bam_left=$2;
mod_bam_right=$3;

# load packages
source load_package.sh -miller -python;

# load configuration variables
source config.sh

# loop over each direction separately and find maximum fork probability value per fork

mod_bam_files=("$mod_bam_left" "$mod_bam_right")
directions=("L" "R")

for cnt in 0 1; do
  # shellcheck disable=SC1010
  # shellcheck disable=SC2016
  mlr --itsv --otsv put '$direction=splitax($detectIndex,"_")[6]' \
    then filter '$direction=="'"${directions[$cnt]}"'"' then cut -f detectIndex\
    then rename detectIndex,alt_read_id "$pause_file" |\
    python get_raw_data_from_modBAM.py --piped-regions --fork-info-in-alt-read-id --get-max "${mod_bam_files[$cnt]}" 
done
