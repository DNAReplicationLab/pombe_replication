#!/usr/bin/env Rscript

# Usage: Rscript plot_correlation.R <input_file> <output_dir>

# Load required package
library(ggplot2)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript plot_correlation.R <input_file> <output_dir>")
}
input_file <- args[1]
output_dir <- args[2]

# Read the combined data file
data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Convert "N/A" to NA and ensure numeric type
data$left_fraction_Puseq[data$left_fraction_Puseq == "N/A"] <- NA
data$left_fraction_Puseq <- as.numeric(data$left_fraction_Puseq)
data$left_fraction <- as.numeric(data$left_fraction)

# Filter out rows with missing values
cor_data <- data[complete.cases(data[, c("left_fraction", "left_fraction_Puseq")]), ]

# Calculate Pearson correlation coefficient
cor_coef <- cor(cor_data$left_fraction, cor_data$left_fraction_Puseq, method = "pearson")
cat("Pearson correlation:", cor_coef, "\n")

# Create the plot
p <- ggplot(cor_data, aes(x = left_fraction, y = left_fraction_Puseq)) +
  geom_point(alpha = 0.4, color = "darkgrey", stroke = 0) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  labs(x = "Pu-seq leftward fork ratios",
       y = "Single molecule\nleftward fork ratios") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14)
  ) +
  annotate("text", 
           x = 0.05, y = 0.95,
           label = paste0("Pearson r = ", round(cor_coef, 3)),
           hjust = 0, size = 5)

# Save the plot
ggsave(filename = file.path(output_dir, "Fig2D", "Fig2D_leftward_fraction_DNAscent_vs_Puseq.svg"),
       plot = p, width = 4.5, height = 2, dpi = 300)
