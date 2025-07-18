require(ggplot2)
require(ggthemes)
require(scales)
options(bitmapType = "cairo")

# goal of program
# ===============
# given many bed files with the CRAC score per interval, arrange intervals in descending order and plot them,
# highlighting the breaks between the files
# Input needed is a directory which contains files named bed_VS_N.bed where N runs from 1 to some number.
# 'VS' means value split, so these are files which have the CRAC score per interval, such that the N = 1
# file has the highest CRAC score intervals and increasing N files have decreasing CRAC score intervals.
# Such a file will be among the outputs of our pause correlation analysis pipeline.
# NOTE: as this script is close to a publication, our focus is on rapidly making the plot and the aesthetics,
#       so some numbers may be hard-coded.

# usage
# =====
# Rscript <script_name>.R /path/to/dir/with/bed/files /path/to/output_file.png

# iterate over the files in the directory and read their contents into a table
dir_path <- commandArgs(trailingOnly = TRUE)[1]
files <- list.files(dir_path, pattern = "bed_VS_[0-9]+.bed", full.names = TRUE)
tbl <- do.call(rbind, lapply(files, function(f) {
  tbl <- read.table(f, header = FALSE, comment.char = "#",
                    col.names = c("contig", "start", "end", "name", "score", "strand", "CRAC_score"))

  # remove any entries in tbl where start and end are equal
  tbl <- tbl[tbl$start != tbl$end,]

  # arrange the CRAC scores in descending order
  tbl <- tbl[order(-tbl$CRAC_score),]

  # add a column for the VS number
  tbl$VS <- as.numeric(gsub("bed_VS_", "", gsub(".bed", "", basename(f))))
  tbl
}))

# create an index column for the plot
tbl$index <- (1:nrow(tbl))

# find the lowest and highest CRAC scores
min_score <- min(tbl$CRAC_score)
max_score <- max(tbl$CRAC_score)
print(log10(min_score))
print(log10(max_score))

# print total number of intervals
print(nrow(tbl))

# identify the boundaries between the files by finding max index for each VS group except the last
tbl$VS_max_index <- ave(tbl$index, tbl$VS, FUN = max)

# get the unique tbl VS max index values
VS_max_index_unique <- unique(tbl$VS_max_index)
VS_max_index_unique_except_last <- VS_max_index_unique[-length(VS_max_index_unique)]

# make the plot
r <- ggplot(tbl, aes(ymax = log10(CRAC_score), ymin = log10(min_score), x = index)) +
  geom_ribbon(fill = "#000000") +
  geom_vline(xintercept = VS_max_index_unique_except_last, colour = "white", linetype = "dashed", size = 4) +
  scale_x_continuous(expand = c(0.012, 0.012), breaks = VS_max_index_unique,
                     guide = guide_axis(angle = 45), limits = c(0, max(VS_max_index_unique))) +
  scale_y_continuous(breaks = c(2,4,6)) +
  theme_classic(base_size = 90) +
  theme(axis.text = element_text(colour = "black"), axis.title.x = element_blank(),
        axis.title.y = element_blank())

ggsave(commandArgs(trailingOnly = TRUE)[2], plot = r, dpi = 400, width = 15, height = 9.1)
