#!/usr/bin/env Rscript

# Goal: Termination count per window and Pu-seq termination frequency, bin the data, and plot the relationship.

library(ggplot2)

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
termination_counts_file <- args[1]
pu_seq_file <- args[2]
output_dir <- args[3]

# Step 1: Load data
file_a <- read.delim(termination_counts_file)      # Termination counts per window
file_b <- read.delim(pu_seq_file)                  # Pu-seq termination frequency

# Step 2: Rename relevant columns for consistency
colnames(file_a)[4] <- "termination_count"
colnames(file_b)[4] <- "termination_frequency"

# Step 3: Merge datasets on genomic coordinates
merged <- merge(
  file_a,
  file_b,
  by = c("contig", "start_window", "end_window"),
  all = TRUE  # Retain unmatched windows
)

# Step 4: Replace missing values with zeros
merged$termination_count[is.na(merged$termination_count)] <- 0
merged$termination_frequency[is.na(merged$termination_frequency)] <- 0

# Step 5: Remove windows where no terminations were observed
merged <- merged[merged$termination_count > 0, ]

# Step 6: Define termination frequency bins
breaks <- c(0, 2, 4, 6, 8, 14)
labels <- c("0-2", "2-4", "4-6", "6-8", "8-14")

# Bin Pu-seq frequency into intervals
merged$freq_bin <- cut(
  merged$termination_frequency,
  breaks = breaks,
  labels = labels,
  include.lowest = FALSE,
  right = FALSE
)

# Step 7: Remove any entries with NA bins (outside bin range)
merged <- merged[!is.na(merged$freq_bin), ]

# Step 8: Plot termination count by frequency bin
plot <- ggplot(merged, aes(x = freq_bin, y = termination_count)) +
  geom_jitter(width = 0.2, alpha = 0.4, color = "darkgrey") +         
  geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA) +        
  theme_classic() +                                                   
  theme(
    text = element_text(size = 14),       # Base text
    axis.text = element_text(size = 14),  # Axis tick labels
    axis.title = element_text(size = 14), # Axis titles
    plot.title = element_text(size = 14), # Title (if used)
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  ) +
  labs(
    x = "Pu-seq termination frequency (%)",
    y = "Single molecule\ntermination count"
  ) +
  ylim(0, 30)

# Step 9: Save the plot to an SVG file
ggsave(
  filename = file.path(output_dir, "Fig2G", "Fig2G_termination_count_vs_frequency.svg"),
  plot = plot,
  width = 4.5,
  height = 2,
  dpi = 300
)
