#!/bin/bash

# the script outputs images and info associated with each image as a pdf file.

# directory which contains package information and configuration information
# i.e. load_package.sh and config.sh
config_dir=..

# load required packages
pwd=$(pwd)
cd $config_dir || exit;
source config.sh
source load_package.sh -miller -latex
cd "$pwd" || exit;

# preamble
# --------
#
# Take a bunch of plots in a folder and put them in a pdf with additional information per plot
#
# sample usage:
# bash <program_name.sh> fit_file png_dir fork_sense_dir additional_png_dir additional_png_prefix feature_prefix
# * only png_dir is required.
# * fork_sense_dir, additional_png_dir, additional_png_prefix, feature_prefix are optional
# * fit_file is also optional, but use /dev/null (a standard linux stand-in for empty file) in its place if unused.
# * meanings of parameters are explained below

# * Basically, png_dir is a directory with many images with filenames in the format plot_readID.png
# * fit_file is tab-separated with headers and each row corresponds to one fork, labelled by the detectIndex.
#   detectIndex is any index that contains the readID in it
#   (fit_file can have other formats, just change the --itsv flag in the miller command below suitably).
# * fork_sense_dir (optional) contains files with left forks, right forks, origins, terminations produced by fork sense.
# * Explanation of "additional_..." parameters:
#   Optionally, additional images associated with a read id can be displayed below each plot_readID.png in the pdf.
#   Set the directory where these additional plots are found, and the prefix in these filenames.
#   only files within this directory starting with the prefix, ending in .png, and contain the readID will be used
#   i.e. the format must be prefix*readID*.png.
# * Explanation of the (optional) "feature_prefix" parameter
#   A read id could be associated with fork features like origins, left/right forks, terminations etc.
#   These can be inferred by this script (normal behaviour), or for the sake of decreasing run time,
#   we can reuse pre-existing files that contain this information (optional behaviour).
#   If such files are available, they must be in the png_dir and must have the filename
#   <feature_prefix>_<readID>. If you don't want to use this feature, simply don't specify this parameter.
#   NOTE: In the line above, angle brackets mean make a substitution.
#   For the document to look nice, all plots should be of the same width.

# Declare directories, files
fit_file=${1:- }
# shellcheck disable=SC2164
png_dir=$(cd "$2"; pwd)
fork_sense_dir=

if [[ -d "${3:- }" ]]; then
  fork_sense_dir=$(cd "$3" || exit; pwd)
fi

additional_png_prefix=${5:- }
additional_png_dir=

if [[ -d "${4:- }" ]] && ! [ "$additional_png_prefix" == " " ];
then
  additional_png_dir=$(cd "$4" || exit; pwd);
  additional_plots_set=1;
else
  additional_plots_set=0;
fi

feature_prefix=${6:- }

# concatenate the read ids while iterating through the plots
read_id_list=

# output latex to a compilation document
{
  echo "\documentclass{article}"
  echo "\usepackage{graphicx}"
  echo "\usepackage{hyperref}"
  echo "\usepackage[a4paper, margin=0.5in]{geometry}"
  echo "\title{Plots of analogue data from selected read ids}"
  echo "\\author{${config[name]} \\\\ ${config[email]}}"
  echo "\date{\today}"
  echo "\begin{document}"
  echo "\maketitle"
  echo "\tableofcontents"
  echo "\pagebreak"

  for FILE in "$png_dir"/plot_*.png
  do

    [ -e "$FILE" ] || continue

    filename=$(echo "$FILE" | rev | cut -d "/" -f 1 | rev)
    read_id_with_extension=$(echo "$filename" | cut -d "_" -f 2)
    read_id=$(echo "$read_id_with_extension" | cut -d "." -f 1)

    # concatenate the read ids while iterating through the plots
    read_id_list+="$read_id\n"

    echo "\section{Read $read_id}"
    echo "\subsection{Plots}"
    echo "\begin{figure}[h!]"
    echo "\centering"
    echo "\includegraphics[scale=0.2]{$FILE}"
    echo "\linebreak"

    if [ "$additional_plots_set" == "1" ]
    then
      for ADDFILE in "$additional_png_dir"/"$additional_png_prefix"*"$read_id"*.png
      do

        [ -e "$ADDFILE" ] || continue
        echo "\includegraphics[scale=0.2]{$ADDFILE}"
        echo "\linebreak"

      done
    fi

    echo "\end{figure}"
    echo "\newpage"

    echo "\subsection{Data}"
    echo "\begin{verbatim}"

    if [ "$feature_prefix" == " " ];
    then
      if [ -f "$png_dir"/plot_data_"$read_id" ]
      then
        < "$png_dir"/plot_data_"$read_id" grep '^#' | fold -w 80
      fi
      bash get_feature_information.sh "$read_id" "$fit_file" "$fork_sense_dir"
    else
      cat "$png_dir"/"$feature_prefix"_"$read_id"
    fi

    echo "{"
    # shellcheck disable=SC2016
    mlr --itsv --ojson --skip-comments filter -s readid="$read_id" '$detectIndex =~ @readid' "$fit_file"
    echo "}"

    echo "\end{verbatim}"
    echo "\pagebreak"
  done

  echo "\section{List of read ids in alphabetical order}"
  echo "\begin{verbatim}"
  echo -e "$read_id_list" | sort
  echo "\end{verbatim}"
  echo "\end{document}"
} > "$png_dir"/compilation.tex

# convert to pdf
# sometimes, pdflatex needs to be run twice. so we are running it twice every time.
# this is not a bug.
pdflatex -output-directory="$png_dir" "$png_dir"/compilation.tex
pdflatex -output-directory="$png_dir" "$png_dir"/compilation.tex
