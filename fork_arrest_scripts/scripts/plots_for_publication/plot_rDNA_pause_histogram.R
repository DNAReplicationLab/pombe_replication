#! /usr/bin/env Rscript
# usage: Rscript scriptName.R somefile.txt optional_suffix
# Written by Anna Rogers and Sathish Thiyagarajan with the help of AI tools (ChatGPT, CoPilot, Gemini etc.).
require(ggplot2)
require(ggthemes)
require(dplyr)
require(tidyr)
require(tidyverse)
require(jsonlite)
require(patchwork)
options(bitmapType = "cairo")

# goal of program
# ===============
# plot a histogram of pauses along rDNA from input data
# NOTE: the input data is binned.

# usage
# =====
# Rscript scriptName.R somefile.json output_dir
# somefile.json must contain the fields:
# - fob1_bedgraph_coarse_left: path to bedgraph file for fob1 with coarse bins on left forks
# - fob1_bedgraph_coarse_right: path to bedgraph file for fob1 with coarse bins on right forks
# - fob1_bedgraph_fine_lead: same as above but with fine bins and for leading strand synthesis
# - fob1_bedgraph_fine_lag: same as above but with fine bins and for lagging strand synthesis
# - wt_bedgraph_lead: path to bedgraph file for wt leading strand
# - wt_bedgraph_lag: path to bedgraph file for wt lagging strand
# - wt_bedgraph_left: path to bedgraph file for wt left forks
# - wt_bedgraph_right: path to bedgraph file for wt right forks
# - wt_ylim: upper y-axis limit for wt plots
# - fob1_ylim: upper y-axis limit for fob1 plots
# - present_as_fold_x: default false. if this is set to true, then we redistribute pauses uniformly
#                      across the rDNA, and then compute the ratio of pauses observed in a bin to the number
#                      of pauses expected in that bin if pauses were uniformly distributed.
#                      In the JSON file, this is a boolean value, so you must set it to true or false (without quotes).
# - fold_x_normalize_together: default false. only used if present_as_fold_x is true.
#                              When present_as_fold_x is set, by default, when we plot say left and right pauses in a
#                              plot together, we normalize each category (left and right) separately.
#                              If this is set to true, then we normalize both left and right pauses together,
#                              so that the ratio of left to right pauses is preserved.

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

# Return the sum of all scores in the data table
sum_score <- function (data_table) {

	# Ensure the input is a data frame
	if (!is.data.frame(data_table)) {
		stop("Input must be a data frame.")
	}

	# Check if the table has exactly 4 columns
	if (ncol(data_table) != 4) {
		stop("Input table must have exactly 4 columns.")
	}

	# Assign temporary column names for easier access
	# Assuming the columns are in the order: contig, start, end, score
	colnames(data_table) <- c("contig", "start", "end", "score")

	# check that the score column is numeric
	if (!is.numeric(data_table$score)) {
		stop("The 'score' column must be numeric.")
	}

	# Calculate sum of all scores
	sum_scores <- sum(data_table$score, na.rm = TRUE)

	return(sum_scores)
}

# Normalize the scores in the data table
normalize_score <- function(data_table, sum_scores = NULL) {

	# Ensure the input is a data frame
	if (!is.data.frame(data_table)) {
		stop("Input must be a data frame.")
	}

	# Check if the table has exactly 4 columns
	if (ncol(data_table) != 4) {
		stop("Input table must have exactly 4 columns.")
	}

	# Calculate sum of all scores if not provided
	if (is.null(sum_scores)) {
		sum_scores <- sum_score(data_table)
	}

	# Assign temporary column names for easier access
	# Assuming the columns are in the order: contig, start, end, score
	colnames(data_table) <- c("contig", "start", "end", "score")

	# Calculate (end - start) for each row
	data_table$length <- data_table$end - data_table$start

	# Calculate sum of all (end - start) values
	sum_lengths <- sum(data_table$length, na.rm = TRUE)

	# Calculate pauses_per_bp
	if (sum_lengths == 0) {
		stop("Sum of (end - start) values is zero, cannot compute pauses_per_bp.")
	}
	pauses_per_bp <- sum_scores / sum_lengths

	# Check if pauses_per_bp is zero to avoid division by zero
	if (pauses_per_bp == 0) {
		stop("pauses_per_bp is zero, cannot normalize scores.")
	} else {
		# Divide each score value by (pauses_per_bp multiplied by that row's (end-start))
		data_table$score_normalized <- data_table$score / (pauses_per_bp * data_table$length)
	}

	# Remove the temporary columns used for calculations, and the column names
	data_table$length <- NULL # Remove the temporary 'length' column
	data_table$score <- NULL # Remove the original 'score' column as we only need the normalized scores
	colnames(data_table) <- NULL

	return(data_table)
}

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

	fob1_bedgraph_coarse_left <- empty_if_null(input_data$fob1_bedgraph_coarse_left)
	fob1_bedgraph_coarse_right <- empty_if_null(input_data$fob1_bedgraph_coarse_right)
	fob1_bedgraph_fine_lead <- empty_if_null(input_data$fob1_bedgraph_fine_lead)
	fob1_bedgraph_fine_lag <- empty_if_null(input_data$fob1_bedgraph_fine_lag)
	wt_bedgraph_lead <- empty_if_null(input_data$wt_bedgraph_lead)
	wt_bedgraph_lag <- empty_if_null(input_data$wt_bedgraph_lag)
	wt_bedgraph_left <- empty_if_null(input_data$wt_bedgraph_left)
	wt_bedgraph_right <- empty_if_null(input_data$wt_bedgraph_right)

	fob1_ylim <- default_if_null(input_data$fob1_ylim, 10)
	wt_ylim <- default_if_null(input_data$wt_ylim, 100)

	present_as_fold_x <- default_if_null(input_data$present_as_fold_x, FALSE)
	if(!is.logical(present_as_fold_x)) {
		stop("present_as_fold_x must be set as a boolean value in the JSON file (true or false without quotes)")
	}

	fold_x_normalize_together <- default_if_null(input_data$fold_x_normalize_together, FALSE)
	if(!is.logical(fold_x_normalize_together)) {
		stop("fold_x_normalize_together must be set as a boolean value in the JSON file (true or false without quotes)")
	}

	if (fold_x_normalize_together && !present_as_fold_x) {
		stop("fold_x_normalize_together can only be set to true if present_as_fold_x is also true")
	}
}

# make the output directory if it doesn't exist
# and change to it
if (!dir.exists(output_dir)) {
	dir.create(output_dir)
}
setwd(output_dir)

# set some plot parameters
fine_x_limits <- c(8.2, 9.137)
coarse_x_limits <- c(0, 9.137)

# ensure that the fine x limits are within the coarse x limits
if (fine_x_limits[1] < coarse_x_limits[1] || fine_x_limits[2] > coarse_x_limits[2]) {
	stop("fine_x_limits must be within coarse_x_limits")
}

if(file.exists(fob1_bedgraph_coarse_left) && file.exists(fob1_bedgraph_coarse_right)){
	#load bedgraph
	tbl <- read.table(fob1_bedgraph_coarse_left, header = FALSE, sep = "")
	tbl2 <- read.table(fob1_bedgraph_coarse_right, header = FALSE, sep = "")
	#tbl is left, tbl2 is right

	if (present_as_fold_x) {
		if (fold_x_normalize_together) {
			# normalize the scores in the bedgraph files together
			sum_scores <- sum_score(tbl) + sum_score(tbl2)
			tbl <- normalize_score(tbl, sum_scores)
			tbl2 <- normalize_score(tbl2, sum_scores)
		} else {
			# normalize the scores in the bedgraph files separately
			tbl <- normalize_score(tbl)
			tbl2 <- normalize_score(tbl2)
		}
	}

	# rename columns (left/right are the values to plot)
	colnames(tbl)[1:4] <- c("contig", "start", "end", "left")
	colnames(tbl2)[1:4] <- c("contig", "start", "end", "right")

	# calculate midpoint of each bar
	tbl$mid_point_kb <- (tbl$start + tbl$end) / 2000
	tbl2$mid_point_kb <- (tbl2$start + tbl2$end) / 2000

	# remove the contig, start, and end from the tables for simplicity during the join
	# (they aren't used during plotting anyway)
	tbl <- tbl[c("mid_point_kb", "left")]
	tbl2 <- tbl2[c("mid_point_kb", "right")]

	# make a new, joined table by the mid_point_kb, new column Direction tells us left or right,
	# column value had the original values
	df_joined <- tbl %>%
	  inner_join(tbl2, by = "mid_point_kb") %>%
	  pivot_longer(cols = !mid_point_kb, names_to = "Direction") %>%
	  filter(mid_point_kb >= coarse_x_limits[1] & mid_point_kb <= coarse_x_limits[2])

	# check that the maximum value is smaller than the y-axis limit, otherwise stop
	if (max(df_joined$value) > fob1_ylim) {
		stop("fob1_ylim is too small for the maximum value in the data")
	}

	# make the plot with 10 pause count (fob1 had max of 5.5) for full length rDNA
	p1 <- ggplot() +
	  geom_bar(data = df_joined, mapping = aes(x = mid_point_kb, y = value, fill=Direction), stat="identity",
			   position=position_dodge(), width = 0.03) +
	  xlab("rDNA coordinate (kb)") +
	  ylab("Fold enrichment of pauses") +
	  xlim(coarse_x_limits) +
	  ylim(c(0, fob1_ylim)) +
	  theme_classic(base_size = 22) +
	  scale_fill_manual(values=c("darkolivegreen", "darkorange")) +
	  theme(axis.text.x=element_text(size=18, color = "black"),
			axis.ticks.length.x = unit(0.15, "cm"),
			axis.text.y=element_text(size=18, color = "black"),
			axis.ticks.length.y = unit(0.15, "cm"),
			axis.title.x = element_text(size = 18),
			axis.title.y = element_text(size = 18),
			legend.position = c(0.15,0.8),
			legend.title = element_text(size = 22),
			legend.background = element_rect(fill = "transparent"),
			panel.grid.minor = element_blank(),
			panel.grid.major = element_blank())

	#to export as a high resolution tiff
	ggsave("fob1_full_rDNA_pauses.tiff",plot = p1,  width = 21.5, height = 10.5, units = "cm", dpi = 300)
} else {
	print("fob1 coarse bedgraph files not found, so skipping fob1 coarse plots")
}

if(file.exists(fob1_bedgraph_fine_lead) && file.exists(fob1_bedgraph_fine_lag)){
	#load bedgraph
	tbl <- read.table(fob1_bedgraph_fine_lead, header = FALSE, sep = "")
	tbl2 <- read.table(fob1_bedgraph_fine_lag, header = FALSE, sep = "")

	if (present_as_fold_x) {
		if (fold_x_normalize_together) {
			# normalize the scores in the bedgraph files together
			sum_scores <- sum_score(tbl) + sum_score(tbl2)
			tbl <- normalize_score(tbl, sum_scores)
			tbl2 <- normalize_score(tbl2, sum_scores)
		} else {
			# normalize the scores in the bedgraph files separately
			tbl <- normalize_score(tbl)
			tbl2 <- normalize_score(tbl2)
		}
	}

	# rename columns (left/right are the values to plot)
	colnames(tbl)[1:4] <- c("contig", "start", "end", "leading")
	colnames(tbl2)[1:4] <- c("contig", "start", "end", "lagging")

	# calculate midpoint of each bar
	tbl$mid_point_kb <- (tbl$start + tbl$end) / 2000
	tbl2$mid_point_kb <- (tbl2$start + tbl2$end) / 2000

	# remove the contig, start, and end from the tables for simplicity during the join
	# (they aren't used during plotting anyway)
	tbl <- tbl[c("mid_point_kb", "leading")]
	tbl2 <- tbl2[c("mid_point_kb", "lagging")]

	# make a new, joined table by the mid_point_kb, new column Direction tells us left or right,
	# column value had the original values
	df_joined <- tbl %>%
	  inner_join(tbl2, by = "mid_point_kb") %>%
	  pivot_longer(cols = !mid_point_kb, names_to = "Strand") %>%
	  filter(mid_point_kb >= fine_x_limits[1] & mid_point_kb <= fine_x_limits[2])

	# check that the maximum value is smaller than the y-axis limit, otherwise stop
	if (max(df_joined$value) > fob1_ylim) {
	  stop("fob1_ylim is too small for the maximum value in the data")
	}

	p2 <- ggplot() +
	  geom_bar(data = df_joined, mapping = aes(x = mid_point_kb, y = value, fill=Strand), stat="identity",
			   position=position_dodge(), width = 0.01) +
	  xlab("rDNA coordinate (kb)") +
	  ylab("Fold enrichment of pauses") +
	  xlim(fine_x_limits) +
	  ylim(c(0, fob1_ylim)) +
	  theme_classic(base_size = 22) +
	  scale_fill_manual(values=c("blue", "red")) +
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

	#to export as a high resolution tiff, keeping same height, but adjusting width
	ggsave("fob1_zoom_rDNA_pauses.tiff",plot = p2,  width = 15.5, height = 10.5, units = "cm", dpi = 300)
} else {
	print("fob1 fine bedgraph files not found, so skipping fob1 fine plots")
}

# now for the wt data, which we want to have both left and right data on the same plot
if(file.exists(wt_bedgraph_left) && file.exists(wt_bedgraph_right)){
	#load bedgraph
	tbl <- read.table(wt_bedgraph_left, header = FALSE, sep = "")
	tbl2 <- read.table(wt_bedgraph_right, header = FALSE, sep = "")

	if (present_as_fold_x) {
		if (fold_x_normalize_together) {
			# normalize the scores in the bedgraph files together
			sum_scores <- sum_score(tbl) + sum_score(tbl2)
			tbl <- normalize_score(tbl, sum_scores)
			tbl2 <- normalize_score(tbl2, sum_scores)
		} else {
			# normalize the scores in the bedgraph files separately
			tbl <- normalize_score(tbl)
			tbl2 <- normalize_score(tbl2)
		}
	}

	# rename columns (lead/lag are the values to plot)
	colnames(tbl)[1:4] <- c("contig", "start", "end", "left")
	colnames(tbl2)[1:4] <- c("contig", "start", "end", "right")

	# calculate midpoint of each bar
	tbl$mid_point_kb <- (tbl$start + tbl$end) / 2000
	tbl2$mid_point_kb <- (tbl2$start + tbl2$end) / 2000

	# remove the contig, start, and end from the tables for simplicity during the join
	# (they aren't used during plotting anyway)
	tbl <- tbl[c("mid_point_kb", "left")]
	tbl2 <- tbl2[c("mid_point_kb", "right")]

	# make a new, joined table by the mid_point_kb, new column Strand tells us lead or lag,
	# column value had the original values
	df_joined <- tbl %>%
	  inner_join(tbl2, by = "mid_point_kb") %>%
	  pivot_longer(cols = !mid_point_kb, names_to = "Direction") %>%
	  filter(mid_point_kb >= coarse_x_limits[1] & mid_point_kb <= coarse_x_limits[2])

	# check that the maximum value is smaller than the y-axis limit, otherwise stop
	if (max(df_joined$value) > wt_ylim) {
		stop("wt_ylim is too small for the maximum value in the data")
	}

	# make the plot with 100 pause count (best used for WT)(wt had max of 80)
	g1 <- ggplot() +
	  geom_bar(data = df_joined, mapping = aes(x = mid_point_kb, y = value, fill=Direction), stat="identity",
			   position=position_dodge(), width = 0.01) +
	  xlab("rDNA coordinate (kb)") +
	  ylab("Fold enrichment of pauses") +
	  xlim(coarse_x_limits) +
	  ylim(c(0, wt_ylim)) +
	  theme_classic(base_size = 22) +
	  scale_fill_manual(values=c("darkolivegreen", "darkorange")) +
	  theme(axis.text.x=element_text(size=18, color = "black"),
			axis.ticks.length.x = unit(0.15, "cm"),
			axis.text.y=element_text(size=18, color = "black"),
			axis.ticks.length.y = unit(0.15, "cm"),
			axis.title.x = element_text(size = 18),
			axis.title.y = element_text(size = 18),
			legend.position = c(0.15,0.8),
			legend.title = element_text(size = 22),
			legend.background = element_rect(fill = "transparent"),
			panel.grid.minor = element_blank(),
			panel.grid.major = element_blank())


	#to export as a high resolution tiff
	ggsave("wt_full_rDNA_pauses.tiff",plot = g1,  width = 21.5, height = 10.5, units = "cm", dpi = 300)

	g1_fob1_levels <- g1 + coord_cartesian(ylim = c(0, fob1_ylim)) + theme(legend.position = "none")
	ggsave("wt_full_rDNA_pauses_fob1_levels.tiff", plot = g1_fob1_levels,  width = 21.5, height = 10.5, units = "cm", dpi = 300)

} else {
	print("wt bedgraph files not found, so skipping wt plots")
}

# now for the wt data, which we want to have both lead and lag data on the same plot
if(file.exists(wt_bedgraph_lead) && file.exists(wt_bedgraph_lag)){
	#load bedgraph
	tbl <- read.table(wt_bedgraph_lead, header = FALSE, sep = "")
	tbl2 <- read.table(wt_bedgraph_lag, header = FALSE, sep = "")

	if (present_as_fold_x) {
		if (fold_x_normalize_together) {
			# normalize the scores in the bedgraph files together
			sum_scores <- sum_score(tbl) + sum_score(tbl2)
			tbl <- normalize_score(tbl, sum_scores)
			tbl2 <- normalize_score(tbl2, sum_scores)
		} else {
			# normalize the scores in the bedgraph files separately
			tbl <- normalize_score(tbl)
			tbl2 <- normalize_score(tbl2)
		}
	}

	# rename columns (lead/lag are the values to plot)
	colnames(tbl)[1:4] <- c("contig", "start", "end", "leading")
	colnames(tbl2)[1:4] <- c("contig", "start", "end", "lagging")

	# calculate midpoint of each bar
	tbl$mid_point_kb <- (tbl$start + tbl$end) / 2000
	tbl2$mid_point_kb <- (tbl2$start + tbl2$end) / 2000

	# remove the contig, start, and end from the tables for simplicity during the join
	# (they aren't used during plotting anyway)
	tbl <- tbl[c("mid_point_kb", "leading")]
	tbl2 <- tbl2[c("mid_point_kb", "lagging")]

	# make a new, joined table by the mid_point_kb, new column Strand tells us lead or lag,
	# column value had the original values
	df_joined <- tbl %>%
	  inner_join(tbl2, by = "mid_point_kb") %>%
	  pivot_longer(cols = !mid_point_kb, names_to = "Strand") %>%
	  filter(mid_point_kb >= fine_x_limits[1] & mid_point_kb <= fine_x_limits[2])

	if (max(df_joined$value) > wt_ylim) {
		stop("wt_ylim is too small for the maximum value in the data")
	}

	#now for a zoomed in version, include origin to end, start at 8.2, end at 9.137 on x-axis
	# make the plot with 100 pause count (best used for WT)(wt had max of 80)
	g2 <- ggplot() +
	  geom_bar(data = df_joined, mapping = aes(x = mid_point_kb, y = value, fill=Strand), stat="identity",
			   position=position_dodge(), width = 0.01) +
	  xlab("rDNA coordinate (kb)") +
	  ylab("Fold enrichment of pauses") +
	  xlim(fine_x_limits) +
	  ylim(c(0, wt_ylim)) +
	  theme_classic(base_size = 22) +
	  scale_fill_manual(values=c("blue", "red")) +
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
	ggsave("wt_zoom_rDNA_pauses.tiff",plot = g2,  width = 15.5, height = 10.5, units = "cm", dpi = 300)
} else {
	print("wt bedgraph files not found, so skipping wt plots")
}

# if all four plots exist, combine them into one
if (exists("p1") && exists("p2") && exists("g1") && exists("g2")) {
	# Hide legend in p2 and p2
	p1 <- p1 + theme(legend.position = "none")
	p2 <- p2 + theme(legend.position = "none")

	g2 <- g2 + annotate("point", x = c(8.855, 8.915, 8.952, 9.027), y = 78, shape = 10)
	g1 <- g1 + annotate("point", x = c(6.747, 7.359, 7.991, 8.111, 8.855, 9.027), y = 78, shape = 10)

	p3 <- g2 + plot_spacer() + g1 + p2 + plot_spacer() + p1 + plot_layout(ncol = 3, widths = c(1, 0.15, 1.3))
	ggsave("rDNA_pauses_combined.tiff", plot = p3, width = 35, height = 21, units = "cm", dpi = 300)
} else {
	print("Not all plots were created, so not combining them")
}
