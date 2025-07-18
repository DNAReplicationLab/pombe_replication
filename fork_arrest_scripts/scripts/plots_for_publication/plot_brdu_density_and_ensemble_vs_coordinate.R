#load packages
library(dplyr)
library(ggplot2)
library(tidyr)
options(bitmapType = "cairo")

# goal of program
# ===============
# plot mean (windowed) BrdU density and one other ensemble measurement against reference coordinate
# The calculation must have already been done; this is just a plotting script.
# NOTE: as this script is close to a publication, our focus is on rapidly making the plot and the aesthetics,
#       so some numbers may be hard-coded. If you are not close to publication and are doing some analysis,
#       then there is likely some other script in the repo which does the same thing as this script but more generally
#       and without refined aesthetics.

# usage
# =====
# < /path/to/brdu/data Rscript <script_name>.R /path/to/ensemble/measurement region column_number \
#     /path/to/output_file.png

# inputs
# ======
# 1. Required BrdU input file must be tab-separated with only two columns and no header.
#    The first column must be the coordinate information in bp and the second column must be the signal information.
#    The region specified on the command line will be used to subset the data (see point 2 below as well).
# 2. Required ensemble measurement input file must be in the bed format.
#    We will only extract data from the region specified on the command line and use data from the specified
#    column number. Region must be specified in the format "chr:start-end" or just "chr".
#    We assume that the BrdU input file contains mean density only on the contig specified in the region,
#    and we only filter the BrdU data based on the start and end of the region.

# get input arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop(paste0("Usage: < /path/to/brdu/data Rscript <script_name>.R /path/to/ensemble/measurement region column_number ",
       "/path/to/output_file.png"))
}

# receive input arguments
ensemble_data_file <- args[1]
region <- args[2]
column_number <- as.numeric(args[3])
output_file <- args[4]

# make an output directory if it doesn't exist
dir.create(dirname(output_file), showWarnings = FALSE)

# read in the data
data <- read.table(file("stdin"), header = FALSE, sep = "\t", comment.char = "#")

# rename columns
colnames(data) <- c("coordinate", "density")

# get the ensemble data
ensemble_data <- read.table(ensemble_data_file, header = FALSE, sep = "\t", comment.char = "#")

# subset the ensemble data
if(grepl(":", region) && grepl("-", region)) {
  contig <- strsplit(region, ":")[[1]][1]
  start <- as.numeric(strsplit(strsplit(region, ":")[[1]][2], "-")[[1]][1])
  end <- as.numeric(strsplit(strsplit(region, ":")[[1]][2], "-")[[1]][2])
} else {
  contig <- region
  start <- 0
  end <- 1e13 # some large number
}
ensemble_data <- ensemble_data %>% filter(V1 == contig & V2 >= start & V3 <= end)

# create a mid-point column and retain only the mid-point and the column of interest
col_of_interest <- paste0("V", column_number)
ensemble_data <- ensemble_data %>% mutate(mid_point = (V2 + V3)/2000) %>% select(mid_point, col_of_interest)

# subset the data as well based on the region start and end
data <- data %>% filter(coordinate >= start & coordinate <= end)

# plot the BrdU density data
p1 <- ggplot(data) + geom_line(aes(x = coordinate/1000, y = density), size = 2, color = "red")

p1 <- p1  +
      labs(x = "Reference coordinate (kb)", y = "BrdU density") +
      theme_classic(base_size = 33) +
      theme(axis.text = element_text(colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_text(colour = "red"),
            axis.ticks.y = element_line(colour = "red"))

# transform the ensemble data to prepare for secondary y-axis
ensemble_data[[col_of_interest]] <- (ensemble_data[[col_of_interest]] - 25) / 35 * 0.75

# Add the secondary y-axis (Trep)
p2 <- p1 +
  geom_line(data = ensemble_data, aes_string(x = "mid_point", y = col_of_interest),
            size = 2, colour = "black") +
  scale_y_continuous(
    name = "BrdU density",
    sec.axis = sec_axis(~ . * 35 / 0.75 + 25, name = "Trep"),
    limits = c(0, 0.75)
  ) +
  theme(
    axis.title.y.right = element_blank(),
    axis.text.y.right = element_text(colour = "black"),
    axis.ticks.y.right = element_line(colour = "black")
  )

ggsave(output_file, plot = p2, dpi = 600, width = 10, height = 4)