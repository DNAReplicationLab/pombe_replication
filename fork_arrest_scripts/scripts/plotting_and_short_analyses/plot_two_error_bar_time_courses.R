require(ggplot2)
require(ggthemes)
require(dplyr)
options(bitmapType = "cairo")

# goal of program
# ===============
# Given two sets of x, y data, calculate error bars for both and plot the error bars only.

# typical usage
# ==============
# Rscript <script_name>.R file_1 file_2 plot.png xcol ycol label_1 label_2 xlabel ylabel xlim error_bar_bin_size
# file_1: first input text file with x and y data, must have columns names
# file_2: second input text file with x and y data, must have columns names
# plot.png: output plot name
# xcol: column name of x data
# ycol: column name of y data
# label_1: label for the first data set
# label_2: label for the second data set
# xlabel: label of x axis
# ylabel: label of y axis
# xlim: of the format num1,num2. These are the lower and upper limits of the xaxis, used for plotting and binning.
# error_bar_bin_size: size of the x axis bin for calculating mean and s.d. for the error bars.

# sample program execution line
# =============================
# let's say we have done measurements on two kinds of cats (persian and burmese),
# and want to plot cat height on y axis and cat weight on x axis, binning cat weight in units of 2 lbs.
# Rscript <script_name>.R /path/to/cat_1.txt /path/to/cat_2.txt height_vs_weight.png weight height persian burmese\
#   "Cat weight (lb)" "Cat height (in)" 2
# NOTE: \ is used to break the line for readability, it is not needed if the line is not broken

# get input arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 11) {
    stop(cat("Usage: Rscript <scriptName.R> <path_to_file_1> <path_to_file_2> <output_plot_file> ",
             "<x_col> <y_col> <label_1> <label_2> <xlabel> <ylabel> <xlimits> <error_bar_bin_size>"))
}

# store file names
file_1 <- args[1]
file_2 <- args[2]

# get data
tbl_1 <- read.table(file_1, header = TRUE, comment.char = "#")
tbl_2 <- read.table(file_2, header = TRUE, comment.char = "#")

# set x and y column names
xcol <- args[4]
ycol <- args[5]

# set labels
label_1 <- args[6]
label_2 <- args[7]

tbl_1$label <- label_1
tbl_2$label <- label_2

# set error bar bin size
error_bar_bin_size <- as.double(args[11])

# combine the two tables
tbl <- rbind(tbl_1, tbl_2)

# set x limits
xlims <- unlist(strsplit(args[10],","))

# make error bars
# create a factor column called x_bin that has levels corresponding to user-set limits and bin size
tbl$x_bin <- as.factor(cut(tbl[,xcol], breaks = seq(as.double(xlims[1]), as.double(xlims[2]), error_bar_bin_size),
                            include.lowest = TRUE, right = FALSE))

# create a new data frame that reports the mean, standard deviation in y, and mean in x grouped by x_bin
# before that, drop any rows that have NA in the x_bin column i.e. data not in the range specified is thrown out
tbl_summary <- tbl %>% filter(!is.na(x_bin)) %>%
                    group_by(x_bin, label) %>%
                    summarise(mean = mean(.data[[ycol]], na.rm = TRUE),
                              sd = sd(.data[[ycol]], na.rm = TRUE),
                              x_mean = mean(.data[[xcol]], na.rm = TRUE))

# create an upper and lower column for the error bars
tbl_summary$upper <- tbl_summary$mean + tbl_summary$sd
tbl_summary$lower <- tbl_summary$mean - tbl_summary$sd

# start plot making
r <- ggplot()

r <- r + geom_errorbar(data = tbl_summary, mapping=aes(x = x_mean, ymin = upper, ymax = lower, color = label)) +
         geom_point(data = tbl_summary, mapping=aes(x = x_mean, y = mean, color = label))

r <- r  +
      xlim(as.double(xlims[1]), as.double(xlims[2])) +
      labs(x = args[8], y = args[9]) +
      theme_bw(base_size = 22)

ggsave(args[3], plot = r, dpi = 400, width = 10, height = 10)