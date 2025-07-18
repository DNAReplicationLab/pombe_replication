require(ggplot2)
require(ggthemes)
require(scales)
options(bitmapType = "cairo")

# goal of program
# ===============
# plot c.d.f of distances between rDNA pauses on the same read
# Input needed is a tab-separated file with columns bin_lo, bin_hi, dist_count
# Such a file will be produced by rDNA_multiple_pause_read_stats.sh
# NOTE: as this script is close to a publication, our focus is on rapidly making the plot and the aesthetics,
#       so some numbers may be hard-coded e.g. rDNA repeat length, etc.

# usage
# =====
# < some_file.txt Rscript <script_name>.R /path/to/output_file.png

tbl <- read.table(file("stdin"), comment.char = "#", header = TRUE)

args <- commandArgs(trailingOnly = TRUE)

repeat_length <- 9.137 # rDNA repeat length in kb

# calculate the midpoint of the bins, and convert lengths to multiples of repeat length
tbl$x <- (tbl$bin_lo + tbl$bin_hi)/(2 * 1000 * repeat_length)

# convert dist_count into a probability distribution
tbl$dist_count <- tbl$dist_count / sum(tbl$dist_count)

# convert distance between pauses into a cumulative distribution
tbl$dist_count <- cumsum(tbl$dist_count) / sum(tbl$dist_count)

r <- ggplot() +
      geom_line(data = tbl, aes(x = x, y = dist_count), color = "black", size = 2) +
      labs(x = "", y = "c.d.f.") +
      scale_x_continuous(breaks = seq(2, 10, by = 2), limits = c(14/repeat_length, 10)) +
      theme_bw(base_size = 40)

ggsave(args[1], plot = r, dpi = 600, width = 10, height = 10)