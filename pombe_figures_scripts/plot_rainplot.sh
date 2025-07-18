#!/bin/bash

#SBATCH --mem=10G
#SBATCH -c 1
#SBATCH -p ei-short
#SBATCH -J lPlotRead
#SBATCH --mail-type=END,FAIL
#SBATCH --time 0:39:59
#SBATCH --constraint=""

# preamble
# --------

# a sample program execution line and its meaning follows

# > bash plot_rainplot.sh sample.mod.bam readID 300 output_dir
# can use sbatch in place of bash.
# > bash plot_rainplot.sh sample.mod.bam readID 300 output_dir mash
# can also optionally provide a prefix like 'mash' in the example above
# > bash plot_rainplot.sh sample.mod.bam readID 300 output_dir mash annotation_file
# can also optionally provide a file called annotation_file in the line above.
# caution: if annotation_file is provided, then a prefix must be provided.
# The annotation file is space-separated and has three columns (without headers):
# start, end, label. start and end refer to locations on the reference genome,
# and label can be 'origin', 'termination', 'leftFork', 'rightFork', or 'pause'.
# > bash plot_rainplot.sh -t 0.6 -b 100 sample.mod.bam readID 300 output_dir mash annotation_file
# * can also optionally provide a threshold in the line above. The threshold is the probability above which
#   a thymidine is called as BrdU. The default threshold is 0.5 and is used if no threshold is provided.
#   Must be a number between 0 and 1.
# * can also optionally provide a window boundary in the line above. The window boundary is the position
#   on the reference genome at which a window is forced to end (remember that window positions are arbitrary to
#   within a window size i.e. the window boundary need not coincide with the start of the read).
#   The default window boundary is -1 and is used if no window boundary is provided.

# the line above extracts data from the read with read id = readID from sample.mod.bam.
# Data is windowed in 300 thymidines after thresholding
# i.e. calling T as BrdU if probability > 0.5.
# If you do not want to show a windowed curve, pass a window size of 0 in the command line invocation.
# The plot and plot data are sent to output_dir and have names like
# plot_readID.png, plot_data_readID.
# If a suffix is specified, say 'mash', then the filenames are mash_readID.png, mash_data_readID
# Plot has two components: the raw data and the windowed data.

# stop execution if any command fails
set -e

# get the threshold from the options if set
while getopts ":t:b:" opt; do
  case $opt in
    t)
      threshold=$OPTARG
      ;;
    b)
      win_boundary=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      ;;
  esac
done

# shift the options out of the command line arguments
shift $((OPTIND - 1))

# set output directory, making it if it doesn't exist
mkdir -p "$4"
op_dir=$(cd "$4"; pwd)

# load packages
source load_package.sh -R -python -samtools -bedtools

# set filenames
temp_file="$op_dir"/"${5:-plot}"_temp_"$2"
data_file="$op_dir"/"${5:-plot}"_data_"$2"
plot_file="$op_dir"/"${5:-plot}"_"$2".png
mod_bam_file="$op_dir"/"${5:-plot}"_"$2".bam

{
  # get information about the read
  bash get_information_from_read_id.sh -o "$mod_bam_file" "$1" "$2";

  # get raw data
  bedtools bamtobed -i "$mod_bam_file" |\
    awk '{print $0 "\t" $4 "_" $1 "_" $2 "_" $3}' |\
    sed '1i\contig\tstart\tend\tread_id\tignore1\tignore2\talt_read_id' |\
    python get_raw_data_from_modBAM.py --piped-regions --alt-read-id-column "$1"

} > "$temp_file"

{

    # output raw data with associated windows of 1 base each
    < "$temp_file" awk '/^#/ {print} !/^#/ {print $1 " " $2 " " $2+1 " " $3 " rawDetect"}'

    # window data
    if ! [ "$3" -eq 0 ];
    then
      < "$temp_file" \
        sed '1i\detectIndex\tposOnRef\tval' |\
        python get_mean_brdU_window.py --window "$3" --thres "${threshold:-0.5}" \
          --forceWinBoundaryAtPos "${win_boundary:--1}" |\
        awk '{print $1 " " $3 " " $4 " " $2 " winDetect"}';
    fi

} > "$data_file"

# plot data
if [ "${6:-flyingTurtleMoons}" == "flyingTurtleMoons" ];
then
  Rscript ./plot_one_read_w_win_or_model_if_needed_pombe_manuscript.R "$data_file" "$plot_file";
else
  Rscript ./plot_one_read_w_win_or_model_if_needed_pombe_manuscript.R "$data_file" "$plot_file" "$6";
fi

# delete temporary file
rm "$temp_file"
