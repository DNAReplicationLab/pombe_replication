#!/bin/bash

#SBATCH --mem=10G
#SBATCH -c 1
#SBATCH -p ei-medium
#SBATCH -J lPlotRead
#SBATCH --mail-type=END,FAIL
#SBATCH --time 0:39:59
#SBATCH --constraint=""

# preamble
# --------

# sample program execution lines and their meanings follow

# > bash limited_plot_read.sh sample.mod.bam readID 300 output_dir
# plots the modification data in raw and 300T-window-averaged representations. can use sbatch in place of bash.
# > bash limited_plot_read.sh sample.mod.bam readID 300 output_dir mash
# can also optionally provide a prefix like 'mash' in the example above
# > bash limited_plot_read.sh sample.mod.bam readID 300 output_dir mash annotation_file
# can also optionally provide a file called annotation_file in the line above.
# caution: if annotation_file is provided, then a prefix ('mash' in the example above) must be provided.
# The annotation file is space-separated and has three columns (without headers):
# start, end, label. start and end refer to locations on the reference genome,
# and label can be 'origin', 'termination', 'leftFork', 'rightFork', or 'pause'.
# > bash limited_plot_read.sh -t 0.6 -b 100 -m C,h -i sample.mod.bam readID 300 output_dir mash annotation_file
# * can also optionally provide a threshold in the line above. The threshold is the probability above which
#   a thymidine is called as BrdU. The default threshold is 0.5 and is used if no threshold is provided.
#   Must be a number between 0 and 1.
# * can also optionally provide a window boundary in the line above. The window boundary is the position
#   on the reference genome at which a window is forced to end (remember that window positions are arbitrary to
#   within a window size i.e. the window boundary need not coincide with the start of the read).
#   The default window boundary is -1 and is used if no window boundary is provided.
# * can also optionally provide a base and a modification code in the line above. The base is the base
#   at which modification information is available and the modification code is the code for the modification;
#   these are 'C' and 'h' respectively in the example above. The default is 'T' and 'T'.
# * use '-i' (i as in reference-independent) if the read you are interested in is unmapped or you want to plot
#   modification data along the basecalled sequence as opposed to the reference sequence.

# this script is a limited version of plot_read.sh.
# the line above extracts data from the read with read id = readID from sample.mod.bam.
# Data is windowed in 300 thymidines after thresholding
# i.e. calling T as BrdU if probability > 0.5.
# If you do not want to show a windowed curve, pass a window size of 0 in the command line invocation.
# The plot and plot data are sent to output_dir and have names like
# plot_readID.png, plot_data_readID.
# If a suffix is specified, say 'mash', then the filenames are mash_readID.png, mash_data_readID
# Plot has two components: the raw data and the windowed data.

# if you want any other functionality, refer to plot_read.sh.

# stop execution if any command fails
set -e

# get the threshold from the options if set
is_ref_independent=0
base=T
code=T
while getopts ":t:b:m:i" opt; do
  case $opt in
    t)
      threshold=$OPTARG
      ;;
    b)
      win_boundary=$OPTARG
      ;;
    m)
      # split the option into base and modification code
      base=$(echo "$OPTARG" | cut -d ',' -f 1)
      code=$(echo "$OPTARG" | cut -d ',' -f 2)
      ;;
    i)
      is_ref_independent=1
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
pwd=$(pwd)
config_dir=..
cd "$config_dir"
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
  if [ "$is_ref_independent" -eq 0 ]; then
    bedtools bamtobed -i "$mod_bam_file" |\
      awk '{print $0 "\t" $4 "_" $1 "_" $2 "_" $3}' |\
      sed '1i\contig\tstart\tend\tread_id\tignore1\tignore2\talt_read_id' |\
      python get_raw_data_from_modBAM.py --piped-regions --alt-read-id-column "$1" --base "$base" --code "$code"
  else
    echo "# do_not_plot_annotations"
    echo "# reporting modification data using basecalled (a.k.a. forward sequence) coordinates"
    python get_raw_data_from_modBAM.py "$mod_bam_file" "$2" ref-independent 0 1000000000000 \
      False "" False "$base" "$code"
      # NOTE: 1000000000000 is just a large number that we use to ensure that the entire read is plotted
  fi

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
cd "$pwd"
if [ "${6:-flyingTurtleMoons}" == "flyingTurtleMoons" ];
then
  Rscript ./plot_one_read_w_win_or_model_if_needed.R "$data_file" "$plot_file";
else
  Rscript ./plot_one_read_w_win_or_model_if_needed.R "$data_file" "$plot_file" "$6";
fi

# delete temporary file
rm "$temp_file"