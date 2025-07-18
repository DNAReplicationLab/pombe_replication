#! /usr/bin/env Rscript
# usage: Rscript scriptName.R somefile.txt optional_suffix
require(ggplot2)
require(ggthemes)
require(jsonlite)
require(dplyr)
require(tidyr)
require(tidyverse)
options(bitmapType = "cairo")

# goal of program
# ===============
# Given leading and lagging bedgraph files containing pause counts,
# plot cumulative pause count versus coordinate.

# usage
# =====
# Rscript scriptName.R somefile.json output_dir
# somefile.json must contain the fields:
# - wt_bedgraph_lead: path to bedgraph file for wt leading strand
# - wt_bedgraph_lag: path to bedgraph file for wt lagging strand

# use bedgraphs whose contents look like the following excerpt (do not assume this is real data and '#'
# is a comment character)

## values multiplied by 1 million for sake of IGV view
## blah blah
## blah
## blah
## blah
#rDNA_22_repeat 2 15 206.016
#rDNA_22_repeat 15 28 274.65
#rDNA_22_repeat 28 41 206.016

# check that the correct number of arguments were provided
if (length(commandArgs(trailingOnly = TRUE)) < 2) {
	stop("Usage: Rscript scriptName.R input_file output_dir")
}

# extract the input file and output directory from the command line arguments
input_file <- commandArgs(trailingOnly = TRUE)[1]
output_dir <- commandArgs(trailingOnly = TRUE)[2]

# read the input file
if (!file.exists(input_file)) {
	stop("Input file does not exist")
} else {
	input_data <- fromJSON(input_file)

	empty_if_null <- function(x) {if (is.null(x)) {return("")} else {return(x)}}
	default_if_null <- function(x, default) {if (is.null(x)) {return(default)} else {return(x)}}

	wt_bedgraph_lead <- empty_if_null(input_data$wt_bedgraph_lead)
	wt_bedgraph_lag <- empty_if_null(input_data$wt_bedgraph_lag)
}

# make the output directory if it doesn't exist
# and change to it
if (!dir.exists(output_dir)) {
	dir.create(output_dir)
}
setwd(output_dir)

# set some plot parameters
coarse_x_limits <- c(0, 9.137)

# now for the wt data, which we want to have both lead and lag data on the same plot
if(file.exists(wt_bedgraph_lead) && file.exists(wt_bedgraph_lag)){
	#load bedgraph
	tbl <- read.table(wt_bedgraph_lead, header = FALSE, sep = "")
	tbl2 <- read.table(wt_bedgraph_lag, header = FALSE, sep = "")

	# rename columns (lead/lag are the values to plot)
	colnames(tbl)[1:4] <- c("contig", "start", "end", "leading_pdf")
	colnames(tbl2)[1:4] <- c("contig", "start", "end", "lagging_pdf")

	# calculate midpoint of each bar
	tbl$mid_point_kb <- (tbl$start + tbl$end) / 2000
	tbl2$mid_point_kb <- (tbl2$start + tbl2$end) / 2000

	# calculate the cumulative sum of the pause counts
	tbl$leading <- cumsum(tbl$leading_pdf)
	tbl2$lagging <- cumsum(tbl2$lagging_pdf)

	# normalize the cumulative sum to the maximum value
	tbl$leading <- tbl$leading / max(tbl$leading)
	tbl2$lagging <- tbl2$lagging / max(tbl2$lagging)

	# remove the contig, start, and end from the tables for simplicity during the join
	# (they aren't used during plotting anyway)
	tbl <- tbl[c("mid_point_kb", "leading")]
	tbl2 <- tbl2[c("mid_point_kb", "lagging")]

	# make a new, joined table by the mid_point_kb, new column Strand tells us lead or lag,
	# column value had the original values
	df_joined <- tbl %>%
	  inner_join(tbl2, by = "mid_point_kb") %>%
	  pivot_longer(cols = !mid_point_kb, names_to = "Strand") %>%
	  filter(mid_point_kb >= coarse_x_limits[1] & mid_point_kb <= coarse_x_limits[2])

	# make the plot
	g2 <- ggplot() +
	  geom_line(data = df_joined, mapping = aes(x = mid_point_kb, y = value, color=Strand), linewidth=1) +
	  xlab("rDNA coordinate (kb)") +
	  ylab("Cumulative distribution") +
	  xlim(coarse_x_limits) +
	  ylim(c(0, 1)) +
	  theme_classic(base_size = 22) +
	  scale_colour_manual(values=c("blue", "red")) +
	  theme(axis.text.x=element_text(size=18, color = "black"),
			axis.ticks.length.x = unit(0.15, "cm"),
			axis.text.y=element_text(size=18, color = "black"),
			axis.ticks.length.y = unit(0.15, "cm"),
			axis.title.x = element_text(size = 18),
			axis.title.y = element_text(size = 18),
			legend.position = c(0.2,0.8),
			legend.title = element_text(size = 22),
			legend.background = element_rect(fill = "transparent"),
			panel.grid.minor = element_blank(),
			panel.grid.major = element_blank())

	#to export as a high resolution tiff
	ggsave("wt_rDNA_pauses_cumulative.tiff",plot = g2,  width = 15.5, height = 10.5, units = "cm", dpi = 300)
} else {
	print("wt bedgraph files not found, so skipping wt plots")
}