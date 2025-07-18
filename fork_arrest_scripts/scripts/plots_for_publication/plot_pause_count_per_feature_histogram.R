#!/usr/bin/env Rscript
# usage: < some_file.txt Rscript <script_name>.R /path/to/output_file.png

# goal of program
# ===============
# given a table with one row per feature (like tRNA1, tRNA2, etc.) and a column for number of pauses,
# plot a histogram of number of pauses per feature.
# NOTE: as this script is close to a publication, our focus is on rapidly making the plot and the aesthetics,
#       so some numbers may be hard-coded. If you are not close to publication and are doing some analysis,
#       then there is likely some other script in the repo which does the same thing as this script but more generally
#       and without refined aesthetics.

# Input file description
#
# * we expect one column of pause counts with no header or column name.
#   each row is interpreted as one feature and the pause count in that feature.
#   each feature could be individual tRNAs, individual genes etc.
# * comment lines are allowed and should start with a '#'

require(ggplot2)
require(ggthemes)
require(dplyr)
options(bitmapType = "cairo")

# load command line arguments if available
if (length(commandArgs(trailingOnly = TRUE)) > 0) {
  args <- commandArgs(trailingOnly = TRUE)
}

if (length(args) != 1) {
  stop("Usage: < some_file.txt Rscript <script_name>.R /path/to/output_file.png", call. = FALSE)
} else {
  df <- read.table(file("stdin"), comment.char = "#", header = FALSE)

  # get a histogram of the number of pauses per feature
  df$X1_count <- as.numeric(df$V1)

  # find the maximum number of pauses in a feature
  max_pauses <- max(df$X1_count)
  x_axis_break <- case_when(
    max_pauses <= 5 ~ 1,
    max_pauses > 5 && max_pauses <= 8 ~ 2,
    max_pauses > 8 && max_pauses <= 15 ~ 5,
    TRUE ~ 1
  )

  p <- ggplot(df, aes(x=X1_count)) +
    geom_histogram(binwidth=1, fill="#003366") +
    scale_y_continuous(trans = scales::pseudo_log_trans(0.5, 10),
                       breaks=c(1, 10, 100)) +
    scale_x_continuous(breaks=seq(0, max_pauses, x_axis_break)) +
    theme_classic(base_size = 44) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())

  ggsave(args[1], plot = p, dpi = 300, width = 6, height = 10)

  # print some summary statistics
  cat("Mean number of pauses per feature: ", mean(df$X1_count), "\n")
  cat("Total number of pauses per feature: ", sum(df$X1_count), "\n")
  cat("Number of features: ", length(df$X1_count), "\n")
  cat("Number of features with pauses: ", sum(df$X1_count > 0), "\n")
}
