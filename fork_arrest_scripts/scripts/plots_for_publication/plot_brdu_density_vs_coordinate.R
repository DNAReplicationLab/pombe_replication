#load packages
library(dplyr)
library(ggplot2)
library(tidyr)
options(bitmapType = "cairo")

# goal of program
# ===============
# plot mean (windowed) BrdU density distribution against reference coordinate.
# The calculation must have already been done; this is just a plotting script.
# NOTE: as this script is close to a publication, our focus is on rapidly making the plot and the aesthetics,
#       so some numbers may be hard-coded. If you are not close to publication and are doing some analysis,
#       then there is likely some other script in the repo which does the same thing as this script but more generally
#       and without refined aesthetics.

# usage
# =====
# < /path/to/brdu/data Rscript <script_name>.R /path/to/output_file.png

# inputs
# ======
# Required input file must be tab-separated with only two columns and no header.
# The first column must be the coordinate information in bp and the second column must be the signal information.

# get input arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: < /path/to/brdu/data Rscript <script_name>.R /path/to/output_file.png")
}

# get output file
output_file <- args[1]

# make an output directory if it doesn't exist
dir.create(dirname(output_file), showWarnings = FALSE)

# read in the data
data <- read.table(file("stdin"), header = FALSE, sep = "\t", comment.char = "#")

# rename columns
colnames(data) <- c("coordinate", "density")

# start plot making
r <- ggplot(data) + geom_line(aes(x = coordinate/1000, y = density), size = 2, color = "black")

r <- r  +
      xlim(-20,20) +
      ylim(0, 0.75) +
      labs(x = "Reference coordinate (kb)", y = "BrdU density") +
      theme_classic(base_size = 44) +
      theme(axis.text = element_text(colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())

ggsave(output_file, plot = r, dpi = 600, width = 10, height = 8)