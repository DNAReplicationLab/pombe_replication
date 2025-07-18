require(ggplot2)
require(ggthemes)
require(dplyr)
options(bitmapType = "cairo")

# goal of program
# ===============
# Given mean, sd and x label (not numerical label but categorical label), plot data with error bars.

# typical usage
# ==============
# < input_data_file Rscript <script_name>.R plot.png xlabel ylabel background_mean,background_sd
# input_data_file: space/tab-separated; with comments (#) and column names, with at least 3 columns: x_label, mean, sd
#                  x_label must be categorical e.g. red, blue etc. instead of numerical e.g. 1, 2, 3 etc.
# plot.png: output plot path
# xlabel: label of x axis
# ylabel: label of y axis
# background_mean,background_sd: mean and sd of background data (e.g. control data) to be plotted as a horizontal line.
#                                Defaults to unset (i.e. no background data plotted).
#                                Format is number,number e.g. 1.2,0.3

# sample program execution line
# =============================
# let's say we have done measurements of height on many kinds of cats (persian, burmese etc.),
# and want to plot cat height mean +- sd on y axis, and cat type on x axis.
# < /path/to/cat/height.tsv Rscript <script_name>.R /path/to/cat/height.png "Cat type" "Cat height (in)"

# get input arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop(cat("Usage: < input_data_file Rscript <script_name>.R plot.png xlabel ylabel background_mean,background_sd\n",
             "background_mean,background sd is optional. For more details, see the script header\n"))
}

# get data
tbl <- read.table(file("stdin"), comment.char = "#", header = TRUE)

# ensure that the input data has at least 3 columns called x_label, mean and sd
if (!("x_label" %in% colnames(tbl) && "mean" %in% colnames(tbl) && "sd" %in% colnames(tbl))) {
    stop(cat("Error: input data must have at least 3 columns called x_label, mean and sd\n"))
}

# create an upper and lower column for the error bars
tbl$upper <- tbl$mean + tbl$sd
tbl$lower <- tbl$mean - tbl$sd

# ensure that the input order of the divisions is preserved
x_label_list <- unique(tbl$x_label)
tbl <- within(tbl, x_label <- factor(x_label, levels=x_label_list))

# start plot making
r <- ggplot()

r <- r + geom_errorbar(data = tbl, mapping=aes(x = x_label, ymin = lower, ymax = upper), color = "blue") +
         geom_point(data = tbl, mapping=aes(x = x_label, y = mean), color = "blue")

# add a horizontal ribbon to show the background mean and sd if supplied
if (length(args) == 4) {

    background_mean_sd <- strsplit(args[4], ",")[[1]]
    background_mean <- as.numeric(background_mean_sd[1])
    background_sd <- as.numeric(background_mean_sd[2])

    # make a copy of tbl and change all upper and lower to background_mean+-background_sd respectively
    tbl_temp_copy <- tbl
    tbl_temp_copy$upper <- background_mean + background_sd
    tbl_temp_copy$lower <- background_mean - background_sd

    r <- r + geom_ribbon(data = tbl_temp_copy, mapping=aes(x = x_label, ymin = lower, ymax = upper, group = 1),
                         fill = "grey", alpha = 0.3)
}

r <- r  +
      labs(x = args[2], y = args[3]) +
      ylim(min(0, min(tbl$lower)), max(tbl$upper)) +
      scale_x_discrete(guide = guide_axis(angle = 45)) +
      theme_bw(base_size = 22)

ggsave(args[1], plot = r, dpi = 400, width = 10, height = 10)