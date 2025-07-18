#load packages
library(dplyr)
library(ggplot2)
library(tidyr)
options(bitmapType = "cairo")

# goal of program
# ===============
# plot expected and observed pause count against reference coordinate
# The calculation must have already been done; this is just a plotting script.
# NOTE: as this script is close to a publication, our focus is on rapidly making the plot and the aesthetics,
#       so some numbers may be hard-coded. If you are not close to publication and are doing some analysis,
#       then there is likely some other script in the repo which does the same thing as this script but more generally
#       and without refined aesthetics.

# usage
# =====
# < /path/to/expected/pause/count/bedgraph Rscript <script_name>.R /path/to/observed/pause/count/bedgraph region \
#     /path/to/output_file.png

# inputs
# ======
# 1. Expected pause count in the usual bedgraph format.
#    The region specified on the command line will be used to subset the data (see point 2 below as well).
# 2. Observed pause count in the usual bedgraph format
#    We will only extract data from the region specified on the command line and use data from the specified
#    column number. Region must be specified in the format "chr:start-end" or just "chr".
#    We assume that the BrdU input file contains mean density only on the contig specified in the region,
#    and we only filter the BrdU data based on the start and end of the region.

# get input arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(paste0("Usage: < /path/to/brdu/data Rscript <script_name>.R /path/to/ensemble/measurement region ",
       "/path/to/output_file.png"))
}

# receive input arguments
observed_data_file <- args[1]
region <- args[2]
output_file <- args[3]

# make an output directory if it doesn't exist
dir.create(dirname(output_file), showWarnings = FALSE)

# read in the data and rename the columns
expected_data <- read.table(file("stdin"), header = FALSE, sep = " ", comment.char = "#")
colnames(expected_data) <- c("contig", "start", "end", "expected_count")

# get the observed data and rename the columns
observed_data <- read.table(observed_data_file, header = FALSE, sep = " ", comment.char = "#")
colnames(observed_data) <- c("contig", "start", "end", "observed_count")

# join the data
data <- inner_join(expected_data, observed_data, by = c("contig", "start", "end"))

# subset the expected and observed data based on the region provided
if(grepl(":", region) && grepl("-", region)) {
  contig_roi <- strsplit(region, ":")[[1]][1]
  start_roi <- as.numeric(strsplit(strsplit(region, ":")[[1]][2], "-")[[1]][1])
  end_roi <- as.numeric(strsplit(strsplit(region, ":")[[1]][2], "-")[[1]][2])
} else {
  contig_roi <- region
  start_roi <- 0
  end_roi <- 1e13 # some large number
}
data <- data %>% filter(contig == contig_roi & start >= start_roi & end <= end_roi)

# create a mid-point column (in kb) and retain only the mid-point and the count column
data <- data %>% mutate(mid_point = (start + end) / 2000) %>% select(mid_point, expected_count, observed_count)

# plot the counts
p <- ggplot(data, aes(x = mid_point)) +
      geom_bar(aes(y = observed_count), stat="identity", colour = "#DDCC77", fill = "#DDCC77") +
      geom_line(aes(y = expected_count), colour = "darkgreen", size = 2, alpha=0.7) +
      scale_y_continuous(expand = c(0, 0), breaks = c(0,5,10),
                         limits = c(0, as.integer(max(data$observed_count) * 1.1))) +
      scale_x_continuous(expand = c(0, 0), limits = c(min(data$mid_point), max(data$mid_point))) +
      theme_classic(base_size = 44) +
      theme(axis.text = element_text(colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_text(colour = "black"),
            axis.ticks.y = element_line(colour = "black"))

# save the plot
ggsave(output_file, plot = p, dpi = 400, width = 35, height = 4)
