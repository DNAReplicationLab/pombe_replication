#!/bin/bash

# goal
# -----
# take the output of create_fork_schematic.sh and convert the red-white gradient to the viridis gradient,
# and the blue obstacle to a red obstacle. This is to ensure the schematic is consistent with the rest of the ms
# and is easy to see.

# usage
# -----
# bash recolour_fork_schematic_to_viridis.sh <input_dir> <output_dir>
#  <input_dir> is the directory where the input files are located. These must have been created by
#    create_fork_schematic.sh.
#  <output_dir> is the directory where the output will be saved. If this directory does not exist, it will be created.
#    We will recolour the image of each frame of the animation, and a gif of the animation,
#    so the filenames of the png and the gif in this directory will be the same as in the input directory.

# stop execution if any command fails
set -e

# check that there are two input arguments
if [ "$#" -ne 2 ]; then
  >&2 echo "Usage: bash recolour_fork_schematic_to_viridis.sh <input_dir> <output_dir>"
  >&2 echo "  <input_dir> is the directory where the input files are located."
  >&2 echo "  <output_dir> is the directory where the output will be saved."
  exit 1;
fi

# load packages
pwd=$(pwd)
cd ..
source load_package.sh -python -imagemagick
cd "$pwd" || exit

# create the schematics
ip_dir="${1:-}"
op_dir="$2"
mkdir -p "$op_dir"

# check that the input directory exists
if [ ! -d "$ip_dir" ]; then
  >&2 echo "Input directory $ip_dir does not exist."
  exit 1;
fi

# recolour the png files corresponding to the frames of the animation
for png_file in "$ip_dir"/*.png; do
  python recolour_fork_schematic_to_viridis.py "$png_file" "$op_dir"/"$(basename "$png_file")"
done

# convert frames to gif
convert -delay 7 -loop 0 "$op_dir"/*_normal*.png -scale 480x480 "$op_dir"/animation_normal_final.gif &
convert -delay 7 -loop 0 "$op_dir"/*_obstacle*.png -scale 480x480 "$op_dir"/animation_obstacle_final.gif &

wait