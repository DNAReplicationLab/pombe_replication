#!/usr/bin/env Rscript
# usage Rscript ./program_name.R input_plot_data output_plot.png

# Input file description
#
# * input_plot_data must have columns: bin_lo,bin_hi,theoretical_count,theoretical_sd,uniform_theoretical_count,
#                                      uniform_theoretical_sd,1_count
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

  p <- ggplot(df, aes(x=bin_lo, y=X1_count)) +
    geom_bar(stat="identity", color="#DDCC77", fill="#DDCC77", position=position_dodge()) +
    geom_ribbon(aes(ymin=theoretical_count - theoretical_sd, ymax=theoretical_count + theoretical_sd),
                color="#117733", fill="#117733", alpha=0.5) +
    geom_ribbon(aes(ymin=uniform_theoretical_count - uniform_theoretical_sd,
                    ymax=uniform_theoretical_count + uniform_theoretical_sd),
                color="#000000", fill="#000000", alpha=0.5) +
    geom_text(aes(label = ifelse(significant == "yes", "*", "")),
      position = position_dodge(width = .9), vjust = -.1, size = 20 / .pt) +
    labs(x="Nearest-neighbour distance (b)", y="Count") +
    theme_classic() +
    scale_fill_manual(values= '#999999') +
    theme_bw(base_size = 22)

  ggsave(args[2], plot = p, dpi = 400, width = 10, height = 10)
}