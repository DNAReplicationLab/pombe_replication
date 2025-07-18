require(ggplot2)
require(ggthemes)
require(dplyr)
options(bitmapType = "cairo")

# goal of program
# ===============
# plot a scatter plot from input data.

# typical usage
# ==============
# < some_file.txt Rscript <script_name>.R column_name_x column_name_y plot.png xlabel ylabel xlim ylim alpha [size]\
#   [colour_column] [error_bar_bin_size] [width] [height] [dpi] [slope_and_intercept]
# NOTE: [] means the argument is optional. If you want to specify the argument, remove the brackets and
#       write the argument. If you don't want to specify the argument, just don't write anything.
# NOTE: The way we've specified optional arguments here, you have to specify all the arguments before the optional
#       argument. For example, if you want to specify the colour_column, you have to specify the size as well.
# some_file.txt: a file with columns separated by one or more space or tabs,
#      with headers out of which two are called column_name_x, column_name_y (see below)
# column_name_x, column_name_y: name of the columns with x and y data
# plot.png: output plot name
# xlabel: label of x axis
# ylabel: label of y axis
# xlim: format is num1,num2 where num1 and num2 are the lower and upper limits of the xaxis. num2 can be set to "auto",
#       in which case the upper limit is set using the maximum value in the column. num1 cannot be set to "auto".
# ylim: same as xlim (above) but for the y axis
# alpha: transparency per point, a number b/w 0 and 1.
# size: (optional) size of the points, default = 1.
# colour_column: (optional), set colour using this column in the input, default = 0 and colour all the points blue.
#                this column must be a column of discrete factors i.e. labels and the like.
#                this column cannot be a column of a continuous number.
# error_bar_bin_size: (optional) size of the x axis bin for calculating mean and s.d. for the error bars.
#                     default = 0 and error bars are not plotted.
# width: (optional) width of the plot (I don't know what the units are), default = 10.
# height: (optional) height of the plot (I don't know what the units are), default = 10.
# dpi: (optional) dots per inch of the plot, default = 400.
# slope_and_intercept: (optional) two values separated by a comma, the slope and intercept of a line to be plotted.
#                      Default = NA,NA i.e. no line is plotted.
#                      This is useful for plotting a line of best fit or a regression line given by y = mx + c.
#                      NOTE: we do not calculate the best-fit slope and intercept,
#                            you have to do that yourself and pass it to the program.

# sample program execution line
# =============================
# let's say we have done measurements on cats, and want to plot cat height on y axis and cat weight on x axis.
# < cat_data.txt Rscript <script_name>.R weight height height_vs_weight.png "Cat weight (g)" "Cat height (in)"\
#     0,5000 0,30 0.1

tbl <- read.table(file("stdin"), comment.char = "#", header = TRUE)

args <- commandArgs(trailingOnly = TRUE)

# store column names
xcol <- args[1]
ycol <- args[2]

# store limits of axes
xlims <- unlist(strsplit(args[6],","))
ylims <- unlist(strsplit(args[7],","))

if(xlims[2] == "auto") {
    xlims[2] <- max(tbl[,xcol]) * 1.1
}

if(ylims[2] == "auto") {
    ylims[2] <- max(tbl[,ycol]) * 1.1
}

if(xlims[1] == "auto" || ylims[1] == "auto") {
    stop("Lower limit of x axis and/or y axis cannot be set to auto")
}

# if error bars are requested, set required flags and parameters
if (length(args) >= 11 && as.double(args[11]) > 0) {
    error_bar_bin_size <- as.double(args[11])
    error_bar_flag <- TRUE
} else {
    error_bar_flag <- FALSE
}

# make a colour blind friendly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

if (length(args) < 8) {
    stop(cat("Usage: < <input_file> Rscript <scriptName.R> <column_name_x> <column_name_y> <output_plot> <xlabel> ",
             "<ylabel> <xlim> <ylim> <alpha> <size (optional, default 1)> <colour_column (optional, default blue)> ",
             "<error_bar_bin_size (optional, default 0)> ",
             "<width (optional, default 10)> <height (optional, default 10)> <dpi (optional, default 400)> ",
             "<slope_and_intercept (optional, default NA,NA)>"))
}

# set a default size for the plot
if (length(args) < 9) {
    args[9] <- 1
}

if (length(args) >= 10){
    if(args[10] == 0){
        set_color_flag <- FALSE
    } else {
        set_color_flag <- TRUE
        # convert the colour column to a factor
        tbl[,args[10]] <- as.factor(tbl[,args[10]])
    }
} else {
    set_color_flag <- FALSE
}

# set width and height and dpi of the plot if specified, or set defaults
if (length(args) >= 12 && as.double(args[12]) > 0) {
    width <- as.double(args[12])
} else {
    width <- 10
}

if (length(args) >= 13 && as.double(args[13]) > 0) {
    height <- as.double(args[13])
} else {
    height <- 10
}

if (length(args) >= 14 && as.integer(args[14]) > 0) {
    dpi <- as.integer(args[14])
} else {
    dpi <- 400
}

# set slope and intercept of line if specified, or set defaults
if (length(args) >= 15 && args[15] != "NA,NA" && args[15] != "") {
    slope_and_intercept <- unlist(strsplit(args[15], ","))
    slope <- as.double(slope_and_intercept[1])
    intercept <- as.double(slope_and_intercept[2])
} else {
    slope <- NA
    intercept <- NA
}

# make error bars if requested
if(error_bar_flag == TRUE) {
    # create a factor column called x_bin that has levels corresponding to user-set limits and bin size
    tbl$x_bin <- as.factor(cut(tbl[,xcol], breaks = seq(as.double(xlims[1]), as.double(xlims[2]), error_bar_bin_size),
                                include.lowest = TRUE, right = FALSE))

    # create a new data frame that reports the mean, standard deviation in y, and mean in x grouped by x_bin
    # before that, drop any rows that have NA in the x_bin column i.e. data not in the range specified is thrown out
    tbl_summary <- tbl %>% filter(!is.na(x_bin)) %>%
                        group_by(x_bin) %>%
                        summarise(mean = mean(.data[[ycol]], na.rm = TRUE),
                                  sd = sd(.data[[ycol]], na.rm = TRUE),
                                  x_mean = mean(.data[[xcol]], na.rm = TRUE))

    # create an upper and lower column for the error bars
    tbl_summary$upper <- tbl_summary$mean + tbl_summary$sd
    tbl_summary$lower <- tbl_summary$mean - tbl_summary$sd

}

# start plot making
r <- ggplot()

if (set_color_flag == FALSE) {
    r <- r + geom_point(data = tbl, mapping = aes_string(x = xcol, y = ycol),
                        alpha = as.double(args[8]), color = "blue", size = as.double(args[9]))
} else {
    # set the colour of the points using the colour_column which is in args[10]
    r <- r + geom_point(data = tbl, mapping = aes_string(x = xcol, y = ycol, color = args[10]),
                        alpha = as.double(args[8]), size = as.double(args[9]))
}

if(error_bar_flag == TRUE) {
    r <- r + geom_errorbar(data = tbl_summary, mapping=aes(x = x_mean, ymin = upper, ymax = lower), color="black") +
             geom_point(data = tbl_summary, mapping=aes(x = x_mean, y = mean), color="black")
}

if(!is.na(slope) && !is.na(intercept)) {
    # add a line with slope and intercept to the plot
    r <- r + geom_abline(slope = slope, intercept = intercept, color = "red", size = 0.5)
}

r <- r  +
      xlim(as.double(xlims[1]), as.double(xlims[2])) +
      ylim(as.double(ylims[1]), as.double(ylims[2])) +
      labs(x = args[4], y = args[5]) +
      theme_bw(base_size = 22)

if (set_color_flag == TRUE) {
    r <- r + scale_fill_manual(values=cbPalette)
}

ggsave(args[3], plot = r, dpi = dpi, width = width, height = height)