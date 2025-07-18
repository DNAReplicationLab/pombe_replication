#!/bin/bash

#SBATCH --mem=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J lPlot1Read
#SBATCH --mail-type=NONE
#SBATCH --mail-user=
#SBATCH --time 0:39:59

# introduction
# --------------

# NOTE: I suggest running plot_n_reads.sh with one read id instead of this script.
#       It is easier to keep track of workflows that way as you use only one script
#       if you want to plot many read ids or just one read id.

# Given a read id, make one plot of raw and windowed data of analogue probabilities,
#                  along with model fit information and fork feature annotations
# In this plot, features on a read id such as forks, pauses, origins, terminations
# are represented. I am not going to describe the representation, as it can change
# as the code is updated. But, it'll be obvious looking at the plot.
# One fork feature out of these can be highlighted.
# To use this feature, set the feature start and feature end suitably (see the section called usage).
# For example, let's say fork sense has called a right fork from 10000 to 20000, among other features.
# Then, to highlight this fork, set feature start and feature end to 10000 and 20000 respectively.
# When the current script looks through what fork sense has called, and it sees a feature whose start
# and end matches the feature start and end specified as inputs to this program, then the feature
# will be highlighted on the plot.

# usage
# --------

# the script takes many inputs, so please read through the comments.

# bash <program_name.sh> [-c] $read_id $mod_bam $feature_start $feature_end $forksense_dir $pause_file $fasta_file\
#     $op_dir [$window_size] [$threshold]
# NOTE: The "\" symbol means the command is continued on the next line.
# NOTE: <> means substitute whatever is within the angle brackets suitably.
# NOTE: [] means the argument is optional. If you want to specify the argument, remove the brackets and
#       write the argument. If you don't want to specify the argument, just don't write anything.
#       To specify the argument, you must specify all arguments to the left of it as well.

# Explanation of inputs
# option -c: If included, plot the reference sigmoid for all forks where pauses were not found.
#            This is useful to see if the sigmoid fits the data well.
#            Not done by default.
#            If you want to specify this option, remove the square brackets and write -c.
# read_id: read id of interest
# mod_bam: the .mod.bam file that contains analogue modification probabilities of the reads of interest
# feature_start: start of the feature of interest, set to fictional numbers or zero to avoid (see introduction)
# feature_end: end of the feature of interest, set to fictional numbers or zero to avoid (see introduction)
# forksense_dir: directory with forksense files that contain origin, left fork, right fork, termination information.
# pause_file: file with pause information. This is a tab separated file with headers.
#             the script uses the columns detectIndex, keep*, pauseSite, pauseDuration,
#             paramStrLeft, paramStrRight, start, end.
#             * means a wildcard i.e. any column that starts with 'keep'.
#             only rows with all keep* values = True are considered as pauses.
#             detectIndex, pauseSite, keep* are required columns.
#             pauseDuration, paramStrLeft, paramStrRight, start, end are optional and will lead to missing
#             plot data/plot features if absent.
#             detectIndex is any string that contains the read id, although it is recommended to use the format
#             readid_contig_start_end_orn_dirn_startFork_endFork where orn = fwd/rev, dirn = L/R, and
#             startFork <= endFork irrespective of fork direction.
#             pauseSite is a position where we've found a potential pause.
#             pauseDuration is the duration of the pause in bases, not in time, although read the script to see if time
#             is ok.
#             keep* are boolean values that indicate whether the pause is genuine or not by different criteria.
#             paramStrLeft, paramStrRight are strings that contain the parameters of model curves fitted to the left
#             and right sides of the pause site respectively.
# fasta_file: reference genome
# op_dir: output directory where all plots are sent.
# window_size: (optional) window size in thymidines to perform averaging for the plot. Default is 300.
# threshold: (optional) threshold above which a thymidine is called as BrdU. Default is 0.5.

# stop execution if any command fails
set -e

# parse options
plot_ref_sgm=0
while getopts "c" opt; do
  case $opt in
    c)
      plot_ref_sgm=1
      ;;
    \?)
      >&2 echo "Invalid option: -$OPTARG"
      exit 1
      ;;
  esac
done

# shift the options off the arguments
shift $((OPTIND-1))

if [ "$#" -lt 8 ]
then
  echo "Incorrect number of arguments! "
  exit
fi

# set configuration directory
config_dir=..

# load packages and configuration variables
pwd=$(pwd)
cd "$config_dir"
source load_package.sh -samtools -bedtools -miller;
source config.sh;
cd "$pwd";

# set input directories, files
read_id=$1
mod_bam=$2
feature_start=$3
feature_end=$4
forksense_dir=
if [[ -d "$5" ]]; then
  forksense_dir=$(cd "$5"; pwd)
fi
pause_file=$6
fa=$7

# set output directory, making it if it doesn't exist
mkdir -p "$8"
op_dir=$(cd "$8"; pwd)

# set window size
window_size=${9:-300}

# set threshold
threshold=${10:-0.5}

# get contig, start, end
contig_start_end=$(samtools view -b -e 'qname=="'"$read_id"'"' "$mod_bam" |\
                   bedtools bamtobed -i stdin |\
                   awk '{print $1 " " $2 " " $3}');

# stop execution if the read couldn't be found.
[ "$contig_start_end" == "" ] && exit 0;

# get feature information to temporary file
bash get_feature_information.sh "$read_id" "$pause_file" "$forksense_dir" > "$op_dir"/features_"$read_id"_temp

{
  # add contig, start, end to feature file
  echo "#Coordinate information: ${contig_start_end}";

  # add more information about the read
  cd "$config_dir";
  bash get_information_from_read_id.sh "$mod_bam" "$read_id";
  cd "$pwd";

  # highlight the feature of interest
  # shellcheck disable=SC2016
  mlr --icsv --ocsv --ifs ' ' --ofs ' ' --headerless-csv-output --skip-comments --implicit-csv-header\
    put -e 'if($1 == '"${feature_start}"' && $2 == '"${feature_end}"'){ $b = 1; } else { $b = 0; }'\
    "$op_dir"/features_"$read_id"_temp;

} > "$op_dir"/features_"$read_id"

# delete temporary file
rm "$op_dir"/features_"$read_id"_temp

# choose the genuine pause site of longest duration out of those in the pause file
# shellcheck disable=SC1010
# shellcheck disable=SC2016
pause_site=$(mlr --itsv --ocsv --ofs ' ' --headerless-csv-output --skip-comments\
  filter -s readid="$read_id" '$detectIndex =~ @readid'\
  then put -e '$b=0;for (k, v in $*){if(k =~ "^keep"){if(v == "False"){$b=1;}}}'\
  then filter '$b==0'\
  then sort -nr pauseDuration\
  then cut -f pauseSite\
  then head -n 1\
  "$pause_file")

# check if the pause file has columns called paramStrLeft, paramStrRight, start, end
# first, extract the header
headerStr=$(mlr --itsv --otsv --skip-comments head -n 1 "$pause_file" | head -n 1)

# check if header has the required columns
pause_site_plot=1
if [[ "$headerStr" =~ "paramStrLeft" ]] && [[ "$headerStr" =~ "paramStrRight" ]] &&\
   [[ "$headerStr" =~ "start" ]] && [[ "$headerStr" =~ "end" ]]
then
  pause_site_plot=1
else
  pause_site_plot=0
fi

# if the option to plot the reference sigmoid is set, then gather the associated parameters
if [ "$plot_ref_sgm" == 1 ]
then

  if [[ "$headerStr" =~ "paramStrNoPause" ]] &&\
   [[ "$headerStr" =~ "start" ]] && [[ "$headerStr" =~ "end" ]] && [[ -f "$pause_file" ]] && [[ -f "$fa" ]]
  then

    # shellcheck disable=SC1010
    # shellcheck disable=SC2016
    param_strings_no_pause=$(mlr --itsv --ocsv --ofs ' ' --headerless-csv-output --skip-comments\
      filter -s readid="$read_id" '$detectIndex =~ @readid'\
      then put -e '$b=0;for (k, v in $*){if(k =~ "^keep"){if(v == "False"){$b=1;}}}'\
      then filter '$b==1'\
      then put '$paramStrNormBdry=$paramStrNoPause . "_" . $start . "_" . $end;
                $paramStrNABdry=($end + 1) . "_NA" '\
      then cut -r -f "paramStr(Norm|NA)Bdry"\
      "$pause_file" | awk 'BEGIN{ORS=" "}{print $0}')

  else

    echo "[ERR] To plot the reference sigmoid, we need the appropriate columns in the pause file and a fasta file." >&2
    exit 1;

  fi

fi

# get parameter strings if need be and plot the read
if [ ! "$pause_site" == "" ] && [ "$pause_site_plot" == 1 ]
then

  # get param strings
  # shellcheck disable=SC1010
  # shellcheck disable=SC2016
  param_strings=$(mlr --itsv --ocsv --ofs ' ' --headerless-csv-output --skip-comments\
    filter -s readid="$read_id" '$detectIndex =~ @readid'\
    then put -e '$b=0;for (k, v in $*){if(k =~ "^keep"){if(v == "False"){$b=1;}}}'\
    then filter '$b==0'\
    then put '$paramStrLeftBdry=$paramStrLeft . "_" . $start . "_" . $pauseSite;
              $paramStrRightBdry=$paramStrRight . "_" . $pauseSite . "_" . $end;
              $paramStrNABdry=($end + 1) . "_NA" '\
    then cut -r -f "paramStr(Left|Right|NA)Bdry"\
    "$pause_file" | awk 'BEGIN{ORS=" "}{print $0}')

  # shellcheck disable=SC2086
  bash plot_read.sh "$mod_bam" "$read_id" $contig_start_end "$op_dir"\
                "$window_size" "$threshold" $pause_site "$fa" $param_strings ${param_strings_no_pause:-} \
                "$op_dir"/features_"$read_id"

else

  if [ "$plot_ref_sgm" == 0 ]
  then
    # shellcheck disable=SC2086
    bash plot_read.sh "$mod_bam" "$read_id" $contig_start_end "$op_dir"\
                  "$window_size" "$threshold" 0 "$op_dir"/features_"$read_id"
  else
    # shellcheck disable=SC2086
    bash plot_read.sh "$mod_bam" "$read_id" $contig_start_end "$op_dir"\
                  "$window_size" "$threshold" 0 "$fa" ${param_strings_no_pause:-} "$op_dir"/features_"$read_id"
  fi

fi
