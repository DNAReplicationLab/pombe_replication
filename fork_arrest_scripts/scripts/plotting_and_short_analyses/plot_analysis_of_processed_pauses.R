#!/usr/bin/env Rscript
require(ggplot2)
require(dplyr)

options(bitmapType = "cairo")

# goal of program
# ===============
# plot several histograms from pause measurements contained in file

# typical usage
# ==============
# Rscript <script_name>.R fit_file.txt output_dir
# fit_file.txt: a file with columns separated by one or more space or tabs,
#      with headers. There are too many required headers to list, so please
#      look through the script.
# output_dir: the directory to which all the plots go

# load command line arguments if available
if (length(commandArgs(trailingOnly = TRUE)) > 0) {
  args <- commandArgs(trailingOnly = TRUE)
}

if (length(args) < 2) {
  stop("Usage: Rscript <script_name>.R fit_file.txt output_dir", call.= FALSE)
}

# load data
fname_fit_results <- args[1]
fit_results_all <- read.table(fname_fit_results, sep = "\t",
                              comment.char = "#", header = TRUE)

# set output directory
op_dir <- args[2]

# set resolution
dpi <- 200

# calculate things like where pauses lie, where forks start etc. using our pause data
if(! "pauseStarts" %in% colnames(fit_results_all))
{
 fit_results_all$pauseStarts <- mapply(function(x, y, w) if (w > 0) x else -y, fit_results_all$leftPieceXmax,
  fit_results_all$rightPieceXmin, fit_results_all$width)
}

if(! "forkStarts" %in% colnames(fit_results_all))
{
  fit_results_all$forkStarts <- mapply(function(x, y, w) if (w > 0) x else -y, fit_results_all$leftPieceXmin,
    fit_results_all$rightPieceXmax, fit_results_all$width)
}

if(! "earlyPieceDurn" %in% colnames(fit_results_all))
{
  fit_results_all$earlyPieceDurn <- mapply(function(x, y, w) if (w > 0) x else y,
    fit_results_all$leftPieceXmax - fit_results_all$leftPieceXmin,
    fit_results_all$rightPieceXmax - fit_results_all$rightPieceXmin,
    fit_results_all$width)
}

if(! "latePieceDurn" %in% colnames(fit_results_all))
{
  fit_results_all$latePieceDurn <- mapply(function(x, y, w) if (w > 0) x else y,
    fit_results_all$rightPieceXmax - fit_results_all$rightPieceXmin,
    fit_results_all$leftPieceXmax - fit_results_all$leftPieceXmin,
    fit_results_all$width)
}

# calculate distance of pause site from the end of the fork it's closest to.
# warning: this quantity will not make any sense for forks where fits haven't worked, so filter suitably later.
if(! "pauseDistFromClosestEnd" %in% colnames(fit_results_all))
{
  fit_results_all$pauseDistFromClosestEnd <- mapply(function(x, y) min(abs(x), abs(y)),
    fit_results_all$pauseSite - fit_results_all$start, fit_results_all$pauseSite - fit_results_all$end
    )
}

if(! "pauseDistFromForkStartNormByForkLen" %in% colnames(fit_results_all))
{
  fit_results_all$pauseDistFromForkStartNormByForkLen <- mapply(function(x, y, w, l) if (w > 0) x/l else y/l,
    fit_results_all$pauseSite - fit_results_all$start, fit_results_all$end - fit_results_all$pauseSite,
     fit_results_all$width, fit_results_all$end - fit_results_all$start)
}

# create a column called all_keep_true which is true if all the columns starting with the word "keep" are set to "True"
fit_results_all$all_keep_true <- apply(select(fit_results_all, starts_with("keep")), 1,
                                       function(x) all(x == "True"))

# plot duration histogram
p <- ggplot(subset(fit_results_all, all_keep_true)) +
      geom_histogram(aes(x = pauseDuration),binwidth = 500, boundary = 0) +
      labs(x = 'Pause duration (b)',y = 'Count') +
      theme_bw(base_size = 22) +
      xlim(0,40000)

# plot pause site histogram
q <- ggplot(subset(fit_results_all, all_keep_true)) +
      geom_histogram(aes(x = pauseStarts),binwidth = 500) +
      labs(x = 'Relative pause location (b)',y = 'Count') +
      theme_bw(base_size = 22)

# plot early-piece duration histogram
r <- ggplot(subset(fit_results_all, all_keep_true)) +
      geom_histogram(aes(x = earlyPieceDurn),binwidth = 200, boundary = 0) +
      labs(x = 'Early piece duration (b)',y = 'Count') +
      theme_bw(base_size = 22)

# plot late-piece duration histogram
s <- ggplot(subset(fit_results_all, all_keep_true)) +
      geom_histogram(aes(x = latePieceDurn),binwidth = 200, boundary = 0) +
      labs(x = 'Late piece duration (b)',y = 'Count') +
      theme_bw(base_size = 22)

# plot fit improvement histogram
t <- ggplot(subset(fit_results_all, all_keep_true)) +
      geom_histogram(aes(x = gofImprovement),binwidth = 0.001, boundary = 0) +
      labs(x = 'Fit improvement',y = 'Count') +
      theme_bw(base_size = 22)

# plot 99 pct pause error calculated using experimental g.o.f. profiles of data with positive duration
v_2 <- ggplot(subset(fit_results_all, pauseDuration > 0 & all_keep_true)) +
      geom_histogram(aes(x = (pauseSite99PctCIHigherUsingExp - pauseSite99PctCILowerUsingExp)/2),
                     binwidth = 100, boundary = 0) +
      labs(x = '99 pct error in pause site (b)',y = 'Count') +
      xlim(0,10000) +
      theme_bw(base_size = 22)

# plot duration histogram for those pauses where the site error is less than 100, calculated using
# expt gof profiles
w_2 <- ggplot(subset(fit_results_all, all_keep_true &
                            (pauseSite99PctCIHigherUsingExp - pauseSite99PctCILowerUsingExp)/2 <= 100)) +
      geom_histogram(aes(x = pauseDuration),binwidth = 500, boundary = 0) +
      labs(x = 'Pause duration, 99 pct site expt err < 100b (b)',y = 'Count') +
      theme_bw(base_size = 22) +
      xlim(0,40000)

# plot start-of-read histogram
a <- ggplot(subset(fit_results_all, all_keep_true)) +
      geom_histogram(aes(x = forkStarts),binwidth = 500) +
      labs(x = 'Relative fork-start location (b)',y = 'Count') +
      theme_bw(base_size = 22)

# plot gof values without cut-and-align
b <- ggplot(subset(fit_results_all, all_keep_true)) +
      geom_histogram(aes(x = gofNoCut),binwidth = 0.005) +
      labs(x = 'Goodness of fit (normalized, no cut)',y = 'Count') +
      theme_bw(base_size = 22)

# plot gof values after cut-and-align
c <- ggplot(subset(fit_results_all, all_keep_true)) +
      geom_histogram(aes(x = gof),binwidth = 0.005) +
      labs(x = 'Goodness of fit (normalized)',y = 'Count') +
      theme_bw(base_size = 22)

# plot gof values without cut-and-align vs pause duration for forks
d <-  ggplot(subset(fit_results_all, all_keep_true), aes(x = gofNoCut, y = pauseDuration)) +
      geom_point(alpha = 0.4, color = "blue") +
      geom_smooth(color = "red", level = 0.99) +
      labs(x = 'Goodness of fit (normalized,no cut)',y = 'Pause duration') +
      theme_bw(base_size = 22)

# plot gof values without cut-and-align vs pause duration for forks
e <-  ggplot(subset(fit_results_all, all_keep_true), aes(x = gofNoCut, y = gofImprovement)) +
      geom_point(alpha = 0.4, color = "blue") +
      geom_smooth(color = "red", level = 0.99) +
      labs(x = 'Goodness of fit (normalized,no cut)',y = 'Fit improvement') +
      theme_bw(base_size = 22)

# plot late piece duration vs gof improvement for forks
f <-  ggplot(subset(fit_results_all, all_keep_true), aes(x = latePieceDurn, y = gofImprovement)) +
      geom_point(alpha = 0.4, color = "blue") +
      geom_smooth(color = "red", level = 0.99) +
      labs(x = 'Late piece duration (b)',y = 'Fit improvement') +
      theme_bw(base_size = 22)

# plot histogram of distances of pause site location from the end of the fork it is closest to
g <- ggplot(subset(fit_results_all, all_keep_true)) +
      geom_histogram(aes(x = pauseDistFromClosestEnd), binwidth = 200, boundary = 0) +
      labs(x = 'Pause site dist from closest fork end (b)',y = 'Count') +
      theme_bw(base_size = 22)

# plot early piece duration vs fork start for forks
h <-  ggplot(subset(fit_results_all, all_keep_true), aes(x = earlyPieceDurn, y = forkStarts)) +
      geom_point(alpha = 0.4, color = "blue") +
      geom_smooth(color = "red", level = 0.99) +
      labs(x = 'Early piece duration (b)',y = 'Relative fork start (b)') +
      theme_bw(base_size = 22)

# plot histogram of maximum fork probability
i <- ggplot(subset(fit_results_all, all_keep_true)) +
      geom_histogram(aes(x = maxForkProb), binwidth = 0.005, boundary = 0) +
      labs(x = 'Maximum fork probability',y = 'Count') +
      theme_bw(base_size = 22)

# plot histogram of distances of pause site location from the start of the fork
j <- ggplot(subset(fit_results_all, all_keep_true)) +
      geom_histogram(aes(x = pauseDistFromForkStartNormByForkLen), binwidth = 0.05, boundary = 0) +
      labs(x = 'Pause dist frm fork start (norm fork len)',y = 'Count') +
      theme_bw(base_size = 22) +
      xlim(0, 1)

# save all graphs
ggsave(
  paste0(op_dir, "/" , "pause_durations.png"),
  plot = p, dpi = dpi)
ggsave(
  paste0(op_dir, "/" , "pause_locations_relative.png"),
  plot = q, dpi = dpi)
ggsave(
  paste0(op_dir, "/" , "early_piece_durations.png"),
  plot = r, dpi = dpi)
ggsave(
  paste0(op_dir, "/" , "late_piece_durations.png"),
  plot = s, dpi = dpi)
ggsave(
  paste0(op_dir, "/" , "fit_improvement.png"),
  plot = t, dpi = dpi)
ggsave(
  paste0(op_dir, "/" , "pause_site_error_99pct_using_exp.png"),
  plot = v_2, dpi = dpi)
ggsave(
  paste0(op_dir, "/" , "pause_duration_w_site_error_99pct_le_100b_using_exp.png"),
  plot = w_2, dpi = dpi)
ggsave(
  paste0(op_dir, "/" , "fork_start_locations_relative.png"),
  plot = a, dpi = dpi)
ggsave(
  paste0(op_dir, "/" , "gof_no_cut_histogram.png"),
  plot = b, dpi = dpi)
ggsave(
  paste0(op_dir, "/" , "gof_histogram.png"),
  plot = c, dpi = dpi)
ggsave(
  paste0(op_dir, "/" , "gof_no_cut_vs_pause_duration.png"),
  plot = d, dpi = dpi)
ggsave(
  paste0(op_dir, "/" , "gof_no_cut_vs_fit_improvement.png"),
  plot = e, dpi = dpi)
ggsave(
  paste0(op_dir, "/" , "late_piece_duration_vs_fit_improvement.png"),
  plot = f, dpi = dpi)
ggsave(
  paste0(op_dir, "/" , "pause_site_min_dist_from_fork_end_histogram.png"),
  plot = g, dpi = dpi)
ggsave(
  paste0(op_dir, "/" , "early_piece_duration_vs_fork_start.png"),
  plot = h, dpi = dpi)
ggsave(
  paste0(op_dir, "/" , "max_fork_prob_histogram.png"),
  plot = i, dpi = dpi)
ggsave(
  paste0(op_dir, "/" , "pause_site_dist_from_fork_start_norm_fork_len_histogram.png"),
  plot = j, dpi = dpi)