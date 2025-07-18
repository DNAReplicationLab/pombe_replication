#!/bin/bash

# the script outputs info associated with the given read id to standard output.

# directory which contains package information and configuration information
# i.e. load_package.sh and config.sh
config_dir=..

# load required packages
pwd=$(pwd)
cd $config_dir || exit;
source config.sh > /dev/null;
source load_package.sh -miller -latex > /dev/null;
cd "$pwd" || exit;

# preamble
# --------
#
# Extract information associated with a read id from forkSense files and fit file.
#
# sample usage:
# bash get_feature_information.sh read_id fit_file fork_sense_dir
#
# meanings of parameters are explained below

# * read_id is the ID of the read of interest
# * fit_file is tab-separated with headers and each row corresponds to one fork, labelled by the detectIndex.
#   detectIndex is a column that contains strings that contain the readID in it
#   (fit_file can have other formats, just change the --itsv flag in the miller command below suitably).
# * fork_sense_dir contains files with left forks, right forks, origins, terminations produced by fork sense.

# get arguments
read_id=${1:- }
fit_file=${2:- }
fork_sense_dir=${3:- }

# check if we have the requisite inputs
{ ! [ "$read_id" == " " ]; } || exit;

# get forksense features
origins="$fork_sense_dir"/origins_DNAscent_forkSense.bed
left_forks="$fork_sense_dir"/leftForks_DNAscent_forkSense.bed
right_forks="$fork_sense_dir"/rightForks_DNAscent_forkSense.bed
terminations="$fork_sense_dir"/terminations_DNAscent_forkSense.bed

# find origins, terminations, and forks corresponding to the read id of interest

# shellcheck disable=SC2016
# shellcheck disable=SC1010
[ -f "$origins" ] && \
  mlr --icsv --ocsv --ifs ' ' --ofs ' ' --implicit-csv-header --skip-comments --headerless-csv-output\
      filter -s readid="$read_id" '$4 =~ @readid' then cut -f 2,3\
      then put '$label="origin"' "$origins"

# shellcheck disable=SC2016
# shellcheck disable=SC1010
[ -f "$left_forks" ]  && \
  mlr --icsv --ocsv --ifs ' ' --ofs ' ' --implicit-csv-header --skip-comments --headerless-csv-output\
      filter -s readid="$read_id" '$4 =~ @readid' then cut -f 2,3\
      then put '$label="leftFork"' "$left_forks"

# shellcheck disable=SC2016
# shellcheck disable=SC1010
[ -f "$right_forks" ] && \
  mlr --icsv --ocsv --ifs ' ' --ofs ' ' --implicit-csv-header --skip-comments --headerless-csv-output\
      filter -s readid="$read_id" '$4 =~ @readid' then cut -f 2,3\
      then put '$label="rightFork"' "$right_forks"

# shellcheck disable=SC2016
# shellcheck disable=SC1010
[ -f "$terminations" ] && \
  mlr --icsv --ocsv --ifs ' ' --ofs ' ' --implicit-csv-header --skip-comments --headerless-csv-output\
      filter -s readid="$read_id" '$4 =~ @readid' then cut -f 2,3\
      then put '$label="termination"' "$terminations"

# find pauses on the read id of interest which are flagged as good pauses in the columns whose name
# starts with "keep"

# shellcheck disable=SC2016
# shellcheck disable=SC1010
[ -f "$fit_file" ] && \
  mlr --itsv --ocsv --ofs ' ' --headerless-csv-output --skip-comments\
      put '$label="pause"'\
      then filter -s readid="$read_id" '$detectIndex =~ @readid'\
      then put -e '$b=0; $c=0; for (k, v in $*){if(k =~ "^keep"){$c+=1; if(v == "False"){$b=1;}}}'\
      then filter '$b==0 && $c >= 0'\
      then cut -f pauseSite,label\
      then put '$pauseSitePlus1=$pauseSite'\
      then cut -o -f pauseSite,pauseSitePlus1,label\
      "$fit_file"

exit 0
