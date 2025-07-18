require(ggplot2)
require(ggthemes)
options(bitmapType = "cairo")

# goal of program
# ===============
# plot a histogram from input data.
# NOTE: the input data is binned.

# typical usage
# ==============
# < some_file.txt Rscript <script_name>.R column_name plot.png xlabel ylabel xlim ylim font_size
# some_file.txt: a file with columns separated by one or more space or tabs,
#      with headers bin_lo, bin_hi, column_name.
# column_name: name of the column with count data
# plot.png: output plot name
# xlabel: label of x axis (optional, default 'x')
# ylabel: label of y axis (optional, default 'Count')
# xlim: (optional) format is num1,num2 where num1 and num2 are the lower and upper limits of the xaxis
# ylim: (optional) format is num1,num2 where num1 and num2 are the lower and upper limits of the yaxis, num2 can be set to "auto"
# font_size: (optional) font size for the plot (default 22)

# histogram x values are the mean of each bin_lo, bin_hi pair
# histogram y values are under column_name

# NOTE: you cannot skip over one optional value to get to another.
# For example, you cannot skip over xlabel, ylabel if you just want to specify xlim.
# You've to specify all three.

tbl <- read.table(file("stdin"), comment.char = "#", header = TRUE)
tbl$x <- (tbl$bin_lo + tbl$bin_hi)/2

args <- commandArgs(trailingOnly = TRUE)

if(is.na(args[3])){ args[3] <- 'x';}
if(is.na(args[4])){ args[4] <- 'Count';}

if(! is.na(args[5])){ xlims <- unlist(strsplit(args[5],",")) }
if(! is.na(args[6])){ ylims <- unlist(strsplit(args[6],",")) }

if(! is.na(args[7])){ font_size <- as.integer(args[7]) } else { font_size <- 22 }

if(exists('ylims') && ylims[2] == "auto") {
    ylims[2] <- max(tbl[,args[1]]) * 1.1
}

r <- ggplot(tbl, aes_string(x = "x", y = args[1])) +
      geom_bar(stat="identity", position=position_dodge(), fill="#003366") +
      {if(exists('xlims')) xlim(as.double(xlims[1]), as.double(xlims[2]))} +
      {if(exists('ylims')) ylim(as.double(ylims[1]), as.double(ylims[2]))} +
      labs(x = args[3], y = args[4]) +
      theme_bw(base_size = font_size)

ggsave(args[2], plot = r, dpi = 400, width = 10, height = 10)
