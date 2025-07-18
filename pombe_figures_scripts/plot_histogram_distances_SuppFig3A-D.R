#!/usr/bin/env Rscript

# Load libraries
library(ggplot2)
library(readr)
library(dplyr)
library(patchwork)

# Read command line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript plot_histogram_distances_SuppFig3A-D.R /path/to/output_dir/Supp_Fig3A-D")
}
output_dir <- args[1]

# Define file info list
file_info <- list(
  data_1 = list(
    path = file.path(output_dir, "closest_distance_to_daigaku_et_al_2015/observed_vs_permutation_summary.csv"),
    title = "Daigaku et al. 2015 (Pu-seq)"
  ),
  data_2 = list(
    path = file.path(output_dir, "closest_distance_to_segurado_et_al_2003/observed_vs_permutation_summary.csv"),
    title = "Segurado et al. 2003"
  ),
  data_3 = list(
    path = file.path(output_dir, "closest_distance_to_heichinger_et_al_2006/observed_vs_permutation_summary.csv"),
    title = "Heichinger et al. 2006"
  ),
  data_4 = list(
    path = file.path(output_dir, "closest_distance_to_hayashi_et_al_2007/observed_vs_permutation_summary.csv"),
    title = "Hayashi et al. 2007"
  )
)

# Plot function
make_plot <- function(file_path, plot_title) {
  df <- read_csv(file_path)

  # Ensure bin is ordered
  df$bin <- factor(df$bin, levels = unique(df$bin), ordered = TRUE)

  ggplot(df, aes(x = bin)) +
    geom_bar(aes(y = obs_count, fill = "Observed"), stat = "identity", show.legend = TRUE) +
    geom_point(aes(y = permutation_mean, color = "Permutation"), size = 0.2, show.legend = TRUE) +
    geom_errorbar(aes(
      ymin = permutation_mean - permutation_sd,
      ymax = permutation_mean + permutation_sd,
      color = "Permutation"
    ), width = 0.1, show.legend = TRUE) +
    labs(title = plot_title, y = "Count", x = "Distance", fill = "", color = "") +
    scale_fill_manual(values = c("Observed" = "grey")) +
    scale_color_manual(values = c("Permutation" = "black")) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 11),
      legend.position = c(0.98, 0.98),
      legend.justification = c("right", "top"),
      legend.background = element_rect(fill = NA, color = NA),
      legend.key = element_rect(fill = NA),
      legend.text = element_text(size = 11)
    ) +
    guides(
      fill = guide_legend(override.aes = list(shape = NA)),
      color = guide_legend(override.aes = list(shape = NA, linetype = 0, size = 4))
    )
}

# Generate all plots
plots <- lapply(file_info, function(info) make_plot(info$path, info$title))

# Combine plots into one figure
combined_plot <- wrap_plots(plots, ncol = 1) +
  plot_annotation(title = "Closest distance from DNAscent initiations to those from other studies")

# Save plot
ggsave(
  filename = file.path(output_dir, "Supp_Fig3A-D_histogram_distances.svg"),
  plot = combined_plot,
  width = 30, height = 30, units = "cm", dpi = 300
)
