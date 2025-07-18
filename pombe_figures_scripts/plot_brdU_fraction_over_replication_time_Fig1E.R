#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[1]

library(dplyr)
library(ggplot2)
library(tidyr)

# Helper function
make_plot <- function(file_names, labels, output_file, facet_labels=NULL) {
  datasets <- lapply(file_names, function(file) {
    df <- read.table(file, header = FALSE, comment.char = "#")
    colnames(df) <- c("contig", "start", "end", "random", "trep", "read_id", "brdU_mean", "thy")
    df <- df %>% filter(thy >= 100)
    grouped <- df %>% group_by(trep) %>% summarise(brdU_means = list(brdU_mean)) %>% unnest(brdU_means)
    intervals <- seq(0, ceiling(max(grouped$trep) / 5) * 5, by = 5)
    grouped$interval <- cut(grouped$trep, breaks = intervals, include.lowest = TRUE)
    return(grouped)
  })

  for (i in seq_along(datasets)) {
    datasets[[i]]$dataset_name <- labels[i]
  }

  combined_data <- do.call(rbind, datasets)
  combined_data$dataset_name <- factor(combined_data$dataset_name, levels = labels)

  plot <- ggplot(combined_data, aes(x = interval, y = brdU_means)) + 
    geom_violin(fill = "darkgray", color = NA, trim = TRUE, scale = "area", bounds = c(0.05, 1)) +
    labs(x = "Replication time (min)", y = "Fraction of BrdU incorporated") +
    theme_classic() +
    facet_wrap(~dataset_name, ncol = 4, labeller = if (!is.null(facet_labels)) labeller(dataset_name = facet_labels) else label_value) +
    ylim(0.05, 1.0) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "gray") +
    stat_summary(fun = median, geom = "point", shape = 18, color = "black", size = 2) +
    scale_x_discrete(breaks = c("(20,25]","(25,30]", "(30,35]", "(35,40]", "(40,45]", "(45,50]", "(50,55]", "(55,60]"),
                     labels = c("20-25", "25-30", "30-35", "35-40", "40-45", "45-50", "50-55", "55-60")) +
    theme(axis.title.x = element_text(size = 13), 
          axis.title.y = element_text(size = 13), 
          axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5), 
          axis.text.y = element_text(size = 11), 
          strip.text = element_text(size = 11))

  ggsave(filename = file.path(output_dir, output_file), plot = plot, width = 10, height = 5, units = "cm", dpi = 300)
}

# --- Top ---
make_plot(
  file_names = file.path(output_dir, c(
    "trep_brdU_mean_per_read_ARY017_10uM-R9-SacCer.txt",
    "trep_brdU_mean_per_read_ARY017_30uM-R9-SacCer.txt",
    "trep_brdU_mean_per_read_ARY017_100uM-R9-SacCer.txt")),
  labels = c("10 µM BrdU", "30 µM BrdU", "100 µM BrdU"),
  output_file = "Fig1E_top.svg"
)

# --- Middle ---
make_plot(
  file_names = file.path(output_dir, c(
    "trep_brdU_mean_per_read_IDS_43-R9-ASM294v2.txt",
    "trep_brdU_mean_per_read_IDS_47-R9-ASM294v2.txt",
    "trep_brdU_mean_per_read_IDS_64-R9-ASM294v2.txt")),
  labels = c("0.5 uM BrdU (IDS_43)", "2 uM BrdU (IDS_47)", "4 uM BrdU (IDS_64)"),
  output_file = "Fig1E_middle.svg",
  facet_labels = c(
    "0.5 uM BrdU (IDS_43)" = "0.5 µM BrdU",
    "2 uM BrdU (IDS_47)" = "2 µM BrdU",
    "4 uM BrdU (IDS_64)" = "4 µM BrdU"
  )
)

# --- Bottom ---
make_plot(
  file_names = file.path(output_dir, c(
    "trep_brdU_mean_per_read_IDS_49-R9-ASM294v2.txt",
    "trep_brdU_mean_per_read_IDS_52-R9-ASM294v2.txt",
    "trep_brdU_mean_per_read_IDS_65-R9-ASM294v2.txt")),
  labels = c("0.5 uM BrdU (IDS_49)", "2 uM BrdU (IDS_52)", "4 uM BrdU (IDS_65)"),
  output_file = "Fig1E_bottom.svg",
  facet_labels = c(
    "0.5 uM BrdU (IDS_49)" = "0.5 µM BrdU",
    "2 uM BrdU (IDS_52)" = "2 µM BrdU",
    "4 uM BrdU (IDS_65)" = "4 µM BrdU"
  )
)
