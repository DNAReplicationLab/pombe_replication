#!/bin/bash

# goal
# -----
# create schematics of moving forks, which we can stitch into a movie

# usage
# -----
# bash create_fork_schematic.sh [-r] <output_dir>
# NOTE: [] denotes optional arguments, so to specify them, please remove the brackets
#   -r: reverse the direction of the gradient
#   <output_dir> is the directory where the output will be saved.
#                The output will be saved as <output_dir>/animation*.png where * is the frame number of the animation

# stop execution if any command fails
set -e

# check that there is at least one input argument
if [ "$#" -lt 1 ]; then
  >&2 echo "Usage: bash create_fork_schematic.sh [-r] <output_dir>"
  >&2 echo "  <output_dir> is the directory where the output will be saved"
  exit 1;
fi

# load packages
pwd=$(pwd)
cd ..
source load_package.sh -python -pdflatex_w_tikz -imagemagick
cd "$pwd" || exit

# set reverse flag
reverse=""
while getopts ":r" opt; do
  case $opt in
    r)
      reverse="-reverse"
      ;;
    \?)
      >&2 echo "Invalid option: $OPTARG"
      >&2 echo "Usage: bash create_fork_schematic.sh [-r] <output_dir>"
      >&2 echo "  <output_dir> is the directory where the output will be saved"
      exit 1;
      ;;
  esac
done

# remove the optional arguments from the input arguments
shift $((OPTIND -1))

# create the schematics
op_dir="$1"
mkdir -p "$op_dir"
tex_file_prefix="$op_dir"/animation
python create_fork_schematic.py $reverse "$tex_file_prefix"_normal
python create_fork_schematic.py $reverse -obstacle "$tex_file_prefix"_obstacle

# convert the tex files corresponding to each frame of the animation to png
for tex_file in "$tex_file_prefix"*.tex; do
  pdf_file=${tex_file%.tex}.pdf
  png_file=${tex_file%.tex}.png
  pdflatex_w_tikz -output-directory="$op_dir" "$tex_file" "$pdf_file"
  convert -density 600 "$pdf_file" "$png_file"
done

# wait for background jobs to finish
wait;

# convert frames to gif
convert -delay 7 -loop 0 "$tex_file_prefix"_normal*.png -scale 480x480 "$tex_file_prefix"_normal_final.gif &
convert -delay 7 -loop 0 "$tex_file_prefix"_obstacle*.png -scale 480x480 "$tex_file_prefix"_obstacle_final.gif &
wait