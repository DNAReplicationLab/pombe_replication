#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Usage: Rscript plot_replication_timing.R <input_bed> <output_pdf> <chromosome> <start:end>")
}

input_bed <- args[1]
output_pdf <- args[2]
chromosome <- args[3]
range <- strsplit(args[4], ":")[[1]]
start_pos <- as.numeric(range[1])
end_pos <- as.numeric(range[2])

library(ggplot2)

data <- read.delim(input_bed, header = FALSE)
filtered_data <- data[data$V1 == chromosome & data$V2 >= start_pos & data$V3 <= end_pos, ]

pdf(output_pdf, width = 6, height = 4)
ggplot(filtered_data, aes(x = V2, y = V5)) +
  geom_line(size = 1) +
  scale_y_reverse() +
  labs(y = "Time (min)", x = paste("Position in", chromosome)) +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 12)
  )
dev.off()
