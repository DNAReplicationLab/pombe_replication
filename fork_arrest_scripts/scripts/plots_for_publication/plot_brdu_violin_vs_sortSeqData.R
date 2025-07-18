#load packages
library(dplyr)
library(ggplot2)
library(tidyr)
options(bitmapType = "cairo")

# Written by Anna Rogers, Isabel Díez Santos, and Sathish Thiyagarajan

# goal of program
# ===============
# plot mean BrdU density distribution vs some externally measured population-level quantity
# The calculation must have already been done; this is just a plotting script.
# Input needed is a tab-separated file with columns
# "contig", "start", "end", "random", "sortSeqData", "read_id", "brdU_mean", "thy" (but no column names in the file)
# NOTE: as this script is close to a publication, our focus is on rapidly making the plot and the aesthetics,
#       so some numbers may be hard-coded. If you are not close to publication and are doing some analysis,
#       then there is likely some other script in the repo which does the same thing as this script but more generally
#       and without refined aesthetics.

# usage
# =====
# < some_file.txt Rscript <script_name>.R /path/to/output_file.png

#read the input file
data <- read.table(file("stdin"), comment.char = "#", header = TRUE)

#parse the output file name
args <- commandArgs(trailingOnly = TRUE)

#renaming columns for ease
colnames(data) <- c("contig", "start", "end", "random", "sortSeq", "read_id", "brdU_mean", "thy")

#remove lines with less than 100 thymidines
filtered_data <- subset(data, thy >= 100)

# Remove rows with NA values in critical columns
filtered_data <- filtered_data %>% drop_na(sortSeq, brdU_mean)

# (can't have copy numbers smaller than 1 or larger than 2, so we just print how many such rows)
print("Number of rows with val < 1")
print(nrow(subset(data, sortSeq < 1)))
print("Number of rows with val > 2")
print(nrow(subset(data, sortSeq > 2)))

# set lines with sortSeq < 1 to 1 and lines with sortSeq > 2 to 2
filtered_data$sortSeq[filtered_data$sortSeq < 1] <- 1
filtered_data$sortSeq[filtered_data$sortSeq > 2] <- 2

# Group the data by sortSeq and collect the values in BrdU_mean for each group
grouped_data <- filtered_data %>% group_by(sortSeq) %>% summarise(brdU_means = list(brdU_mean))

# All of these comments from Isabel’s notes, Convert the tibble to long format using tidyr::unnest().
# ggplot2 expects the data to be in a "long" format, where each row corresponds to an individual observation.
# However, in this case, the data is in a "wide" #format, where some rows have multiple values in a single column.
# To resolve this issue, you need to transform the data into a long format before plotting.
# You can use the tidyr package to achieve this.
grouped_data_long <- grouped_data %>% unnest(brdU_means)

#Define the number of sortSeq groups
intervals <- c(1, 1.2, 1.4, 1.6, 1.8, 2)

#add the internal to the data
grouped_data_long$interval <- cut(grouped_data_long$sortSeq, breaks = intervals, include.lowest = TRUE)

# make the plot
p <- ggplot(grouped_data_long, aes(x = interval, y = brdU_means)) +
  geom_violin(fill = "darkgray", color = "NA", trim = TRUE, scale = "area", bounds=c(0.1, 1)) +
  labs(x = "Sort-seq copy number", y = "BrdU density") +
  theme_classic(base_size = 38)

p <- p + theme(
  axis.title.x = element_text(color = "black", size = 18),
  axis.title.y = element_text(color = "black", size = 18),
  axis.text = element_text(color = "black", size = 16),
  axis.line = element_line(colour = "black", linewidth = 1, linetype = "solid"),
  axis.ticks = element_line(color = "black", linewidth = 0.7)) +
  ylim(0.1, 1) +
  scale_x_discrete(breaks=c("[1,1.2]", "(1.2,1.4]", "(1.4,1.6]", "(1.6,1.8]", "(1.8,2]"),
                   labels=c("1.0-1.2", "1.2-1.4", "1.4-1.6", "1.6-1.8", "1.8-2.0")) +
  geom_hline(yintercept = 0.1, linetype="dashed", color="gray") +
  stat_summary(fun = median, geom = "point", shape=18, color="black", size=4)

ggsave(args[1], plot = p, dpi = 400, width = 6.5, height = 5)
