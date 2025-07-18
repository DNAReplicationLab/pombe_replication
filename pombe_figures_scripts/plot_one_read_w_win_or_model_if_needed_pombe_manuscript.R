#!/usr/bin/env Rscript
# usage: Rscript ./program_name.R plot_data plot.png plot_elements_this_is_optional

# Load required libraries
library(ggplot2)
library(ggthemes)
library(dplyr)
options(bitmapType = "cairo")

# Set dpi and dimensions of output plot
dpi <- 600
width <- 18
height <- 12
base_size <- 60

# Set default colour of the plot of windowed detect data
colour_win_detect <- "black"

# Load command line arguments if available
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript ./<script_name>.R plot_data plot.png plot_elements_this_is_optional", call. = FALSE)
}

# Process the annotation file if it exists
annotations <- NULL
if (length(args) >= 3 && file.exists(args[3]) && file.size(args[3]) > 0) {
  # Read annotations
  annotations <- read.table(args[3], header = FALSE, comment.char = "#")

  if (ncol(annotations) == 3) {
    colnames(annotations) <- c("start", "end", "label")
    annotations$size <- 1
  } else if (ncol(annotations) == 4) {
    colnames(annotations) <- c("start", "end", "label", "size")
  } else {
    stop("plot_elements_this_is_optional must have three or four columns", call. = FALSE)
  }

  # Convert positions to kb
  annotations$start <- annotations$start / 1000
  annotations$end <- annotations$end / 1000

  # Set y-coordinates based on annotation type
  annotations$y_start <- 0.98
  annotations$y_end <- 0.98

  # Assign colors based on label
  annotations$colour <- "black"  
  annotations$colour[annotations$label == "rightFork"] <- "#F8766D"  
  annotations$colour[annotations$label == "leftFork"] <- "#619CFF"   
  annotations$colour[annotations$label == "termination"] <- "grey"
  annotations$colour[annotations$label == "pause"] <- "black"
  annotations$colour[annotations$label == "origin"] <- "black"

  # Adjust y-coordinates for specific annotations
  annotations$y_start[annotations$label %in% c("pause", "termination", "origin")] <- 1.00
  annotations$y_end[annotations$label %in% c("termination", "origin")] <- 1.00
  annotations$y_end[annotations$label == "pause"] <- 0.96
}

# Load data
read_data <- read.table(args[1], header = FALSE, comment.char = "#")
colnames(read_data) <- c("id", "start", "end", "val", "label")

# Check that the id column has only one value
if (length(unique(read_data$id)) > 1) {
  stop("id column has more than one value! Script behaviour is undefined in this case", call. = FALSE)
}

# Convert positions to kb
read_data$start <- read_data$start / 1000
read_data$end <- read_data$end / 1000

# Identify which labels are present
indices <- c()
if ('rawDetect' %in% read_data$label) indices <- c(indices, 1)
if ('winDetect' %in% read_data$label) indices <- c(indices, 2)

# Plot the data
plot1 <- ggplot() +
  geom_point(data = subset(read_data, label == "rawDetect"),
             aes(x = start, y = val), shape = 'circle', alpha = 0.2, colour = "#888888", show.legend = TRUE) +
  geom_segment(data = subset(read_data, label == "winDetect"),
               aes(x = start, y = val, xend = end, yend = val), colour = colour_win_detect, linewidth = 2, show.legend = FALSE) +
  geom_step(data = subset(read_data, label == "winDetect"),
            aes(x = start, y = val), colour = colour_win_detect, linewidth = 2, show.legend = TRUE) +
  xlab("Reference coordinate (kb)") +
  ylab("Probability of BrdU") +
  ylim(c(0, 1)) +
  scale_colour_manual(name = NULL,
                      values = c("rawDetect" = "#888888", "winDetect" = colour_win_detect),
                      limits = c('rawDetect', 'winDetect'),
                      breaks = indices,
                      labels = c("rawDetect" = "Raw data", "winDetect" = "Windowed data"),
                      guide = guide_legend(override.aes = list(linetype = c("blank", "solid")[indices]))) +
  theme_classic(base_size = base_size) +
  theme(legend.position = "none", panel.grid.minor = element_blank(), panel.grid.major = element_blank())

# Add annotations if they exist
if (!is.null(annotations)) {
  plot1 <- plot1 +
    geom_segment(data = subset(annotations, label == "rightFork"),
                 aes(x = start, y = y_start, xend = end, yend = y_end, colour = colour),
                 arrow = arrow(length = unit(0.8, "cm"), type = "closed", angle = 10),
                 linewidth = 4, alpha = 0.5) +
    geom_segment(data = subset(annotations, label == "leftFork"),
                 aes(x = end, y = y_start, xend = start, yend = y_end, colour = colour),
                 arrow = arrow(length = unit(0.8, "cm"), type = "closed", angle = 10),
                 linewidth = 4, alpha = 0.5) +
    geom_segment(data = subset(annotations, label == "pause"),
                 aes(x = start, y = y_start, xend = start, yend = y_end),
                 colour = "black", linewidth = 4) +
    geom_segment(data = subset(annotations, label == "origin"),
                 aes(x = start, y = y_start, xend = end, yend = y_end),
                 colour = "black", linewidth = 4) +
    geom_segment(data = subset(annotations, label == "termination"),
                 aes(x = start, y = y_start, xend = end, yend = y_end),
                 colour = "grey", linewidth = 4) +
    scale_color_identity()
}

# Save the plot
ggsave(args[2], plot = plot1, dpi = dpi, width = width, height = height)
