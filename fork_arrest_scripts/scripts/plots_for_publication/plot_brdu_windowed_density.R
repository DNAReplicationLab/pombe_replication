require(ggplot2)
require(ggthemes)
options(bitmapType = "cairo")

# goal of program
# ===============
# plot windowed brdu density distribution
# Input needed is a space-separated file with columns bin_lo, bin_hi, mean_brdU_count
# Such a file will be among the outputs of dnascent_qc
# NOTE: as this script is close to a publication, our focus is on rapidly making the plot and the aesthetics,
#       so some numbers may be hard-coded. If you are not close to publication and are doing some analysis,
#       then there is likely some other script in the repo which does the same thing as this script but more generally
#       and without refined aesthetics.

# usage
# =====
# < some_file.txt Rscript <script_name>.R /path/to/output_file.png

tbl <- read.table(file("stdin"), comment.char = "#", header = TRUE)
tbl$x <- (tbl$bin_lo + tbl$bin_hi)/2

# normalize the mean_brdU_count column by its sum
tbl$mean_brdU_count <- tbl$mean_brdU_count/sum(tbl$mean_brdU_count)

args <- commandArgs(trailingOnly = TRUE)

r <- ggplot(tbl, aes(x = x, y = mean_brdU_count)) +
      geom_bar(stat="identity", position=position_dodge(), fill="#003366") +
      scale_x_continuous(expand = c(0.0015, 0.0015), limits = c(0, 1.1), breaks = seq(0, 1, 0.25)) +
      scale_y_continuous(expand = c(0.0015, 0.0015), limits = c(0, 0.04), breaks = seq(0, 0.04, 0.01)) +
      xlab("Windowed BrdU density") +
      ylab("Probability density") +
      theme_classic(base_size = 38) +
      theme(axis.text = element_text(colour = "black"))

ggsave(args[1], plot = r, dpi = 400, width = 10, height = 8)
