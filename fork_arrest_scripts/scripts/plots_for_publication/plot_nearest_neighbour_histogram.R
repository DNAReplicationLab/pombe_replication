#!/usr/bin/env Rscript
# usage Rscript ./program_name.R input_plot_data output_plot.png

# goal of program
# ===============
# plot expected and observed nearest-neighbour distribution of pauses along a reference genome
# NOTE: as this script is close to a publication, our focus is on rapidly making the plot and the aesthetics,
#       so some numbers may be hard-coded. If you are not close to publication and are doing some analysis,
#       then there is likely some other script in the repo which does the same thing as this script but more generally
#       and without refined aesthetics.

# Input file description
#
# * input_plot_data must have columns: bin_lo,bin_hi,theoretical_count,uniform_theoretical_count,1_count
# * Column names are required and columns are separated by commas.
# * comment lines begin with #
# * The columns are related to the number of pauses separated by some distance in the genome,
#   and how many pauses are expected from a couple of theoretical models.

require(ggplot2)
require(ggthemes)
options(bitmapType = "cairo")

# load command line arguments if available
if (length(commandArgs(trailingOnly = TRUE)) > 0) {
  args <- commandArgs(trailingOnly = TRUE)
}

if (length(args) < 2) {
  stop("Usage: Rscript ./<script_name>.R input_plot_data output_plot.png", call.= FALSE)
}

# process the annotation file if it exists
if (length(args) == 2 && file.size(args[1]) > 0) {

  df <- read.table(args[1], comment.char = "#", header = TRUE, sep=",")
  df$significant <- ifelse(abs(df$X1_count - df$theoretical_count) >= 3 * df$theoretical_sd, "yes", "no")

  # get midpoint of intervals and convert bp to kbp for plotting
  df$X <- (df$bin_lo + df$bin_hi) / (2*1000)

  p <- ggplot(df, aes(x=X, y=X1_count)) +
    geom_bar(stat="identity", color="#DDCC77", fill="#DDCC77", position=position_dodge()) +
    geom_line(aes(y=theoretical_count), color="#117733", linewidth=2) +
    geom_line(aes(y=uniform_theoretical_count), color="#000000", linewidth=2) +
    scale_fill_manual(values= '#999999') +
    theme_classic(base_size = 44) + theme(axis.title.x = element_blank(), axis.title.y = element_blank())

  ggsave(args[2], plot = p, dpi = 300, width = 10, height = 10)
}
