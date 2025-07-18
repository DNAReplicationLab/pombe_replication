#!/usr/bin/env Rscript

# Usage: Rscript plot_count_vs_efficiency.R <initiations_file> <pu_seq_file> <output_dir>

library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
initiations_file <- args[1]
puseq_file <- args[2]
output_dir <- args[3]

# Load Pu-seq origin data
Puseq_data <- read.delim(puseq_file, header = TRUE)

# Ensure required columns exist
if (!all(c("contig", "Position", "Efficiency..") %in% colnames(Puseq_data))) {
  stop("Expected columns: contig, Position, Efficiency..")
}

Puseq_data <- Puseq_data %>%
  transmute(contig = contig,
            Puseq_origin_position = Position,
            efficiency = `Efficiency..`)

# Function to generate 2000-bp window around origin
generate_windows <- function(origin_position, window_size = 2000) {
  half <- window_size %/% 2
  start <- origin_position - half
  end <- origin_position + half - 1
  tibble(start_window = start, end_window = end)
}

# Apply window generation
Puseq_windows <- Puseq_data %>%
  rowwise() %>%
  mutate(window = list(generate_windows(Puseq_origin_position))) %>%
  unnest(cols = c(window))

# Load initiation midpoints
origins <- read.table(initiations_file, header = FALSE, sep = "\t",
                      col.names = c("contig", "origin_midpoint"))

# Count overlaps
overlap_counts <- mapply(function(start, end, contig) {
  sum(origins$contig == contig & origins$origin_midpoint >= start & origins$origin_midpoint <= end)
}, Puseq_windows$start_window, Puseq_windows$end_window, Puseq_windows$contig)

# Combine counts with metadata
result <- Puseq_windows %>%
  mutate(count = overlap_counts,
         efficiency = Puseq_data$efficiency)

# Bin efficiency
data2 <- result %>%
  mutate(efficiency_group = cut(efficiency,
                                breaks = c(0, 10, 20, 30, 40, 50, 60, 100),
                                labels = c("0-10", "10-20", "20-30", "30-40",
                                           "40-50", "50-60", "60-100"),
                                include.lowest = TRUE)) %>%
  na.omit()

# Plot
plot <- ggplot(data2, aes(x = efficiency_group, y = count)) +
  geom_jitter(alpha = 0.5, size = 1, color = "darkgrey") +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  theme_classic() +
  labs(x = "Pu-seq efficiency (%)", y = "Single molecule\ninitiation count") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  ylim(0, 100)

# Save
ggsave(filename = file.path(output_dir, "Fig2F", "Fig2F_count_vs_efficiency.svg"),
       plot = plot, width = 4.5, height = 2, dpi = 300)
