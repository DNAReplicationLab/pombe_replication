require(ggplot2)
require(ggthemes)
require(dplyr)
options(bitmapType = "cairo")

# goal of program
# ===============
# Given mean and x label (not numerical label but categorical label), plot a scatter plot of the data.

# typical usage
# ==============
# < input_data_file Rscript <script_name>.R plot.png xlabel ylabel
# input_data_file: space/tab-separated; with comments (#) and column names, with at least 2 columns: x_label, mean
#                  x_label must be categorical e.g. red, blue etc. instead of numerical e.g. 1, 2, 3 etc.
# plot.png: output plot path
# xlabel: label of x axis
# ylabel: label of y axis

# sample program execution line
# =============================
# let's say we have done measurements of height on many kinds of cats (persian, burmese etc.),
# and want to plot cat height mean on y axis, and cat type on x axis.
# < /path/to/cat/height.tsv Rscript <script_name>.R /path/to/cat/height.png "Cat type" "Cat height (in)"

# get input arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop(cat("Usage: < input_data_file Rscript <script_name>.R plot.png xlabel ylabel\n"))
}

# get data
tbl <- read.table(file("stdin"), comment.char = "#", header = TRUE)

# ensure that the input data has at least 2 columns called x_label and mean
if (!("x_label" %in% colnames(tbl) && "mean" %in% colnames(tbl))) {
    stop(cat("Error: input data must have at least 2 columns called x_label and mean\n"))
}

# ensure that the input order of the divisions is preserved
x_label_list <- unique(tbl$x_label)
tbl <- within(tbl, x_label <- factor(x_label, levels=x_label_list))

# start plot making
r <- ggplot()

r <- r + geom_point(data = tbl, mapping=aes(x = x_label, y = mean), color = "blue")

r <- r  +
      labs(x = args[2], y = args[3]) +
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      theme_bw(base_size = 22)

ggsave(args[1], plot = r, dpi = 400, width = 10, height = 10)