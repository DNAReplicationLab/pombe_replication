require(ggplot2)
require(ggthemes)
require(scales)
require(ggsankey)
require(dplyr)
options(bitmapType = "cairo")

# goal of program
# ===============
# Given set A and set B each with many bed files with a CRAC score per interval,
# match bed intervals between set A and set B and
# - plot the two CRAC scores of the matched intervals (log-log).
# - optionally make a sankey diagram of the value splits.
# Input needed is two directories which contains files named bed_VS_N.bed where N runs from 1 to some number.
# 'VS' means value split, so these are files which have the CRAC score per interval, such that the N = 1
# file has the highest CRAC score intervals and increasing N files have decreasing CRAC score intervals.
# Such files will be among the outputs of our pause correlation analysis pipeline.
# We will join the two sets of bed files by matching the intervals between the two sets and plot CRAC score 1
# on the X axis and CRAC score 2 on the Y axis (see command below for usage).
# NOTE: as this script is close to a publication, our focus is on rapidly making the plot and the aesthetics,
#       so some numbers may be hard-coded.

# usage
# =====
# Rscript <script_name>.R /path/to/dir_1/with/bed/files /path/to/dir_2/with/bed/files \
#   /path/to/output_file.png [/path/to/output_sankey.png]
# NOTE: \ means the command continues on the next line.
# NOTE: [] means the sankey diagram is optional. To make it, provide the output file path, and remove the [].

read_bed_file <- function (f) {
  tbl <- read.table(f, header = FALSE, comment.char = "#",
                    col.names = c("contig", "start", "end", "name", "score", "strand", "CRAC_score"))

  # remove any entries in tbl where start and end are equal
  tbl <- tbl[tbl$start != tbl$end,]

  # remove any entries in tbl where CRAC_score is NA or zero
  tbl <- tbl[!is.na(tbl$CRAC_score) & tbl$CRAC_score != 0,]

  # append chr, start, end to name
  tbl$name <- paste(tbl$name, tbl$contig, tbl$start, tbl$end, sep = "_")

  # check that the names are unique or raise an error
  if (length(unique(tbl$name)) != nrow(tbl)) {
      stop("The names in the file ", f, " are not unique.")
  }

  # add a column for the VS number
  tbl$VS <- paste0("Q", gsub("bed_VS_", "", gsub(".bed", "", basename(f))))

  tbl
}

# iterate over the files in the directory and read them into a single table for set A
dir_path_a <- commandArgs(trailingOnly = TRUE)[1]
files_a <- list.files(dir_path_a, pattern = "bed_VS_[0-9]+.bed", full.names = TRUE)
tbl_a <- do.call(rbind, lapply(files_a, read_bed_file))

# iterate over the files in the directory and read them into a single table for set B
dir_path_b <- commandArgs(trailingOnly = TRUE)[2]
files_b <- list.files(dir_path_b, pattern = "bed_VS_[0-9]+.bed", full.names = TRUE)
tbl_b <- do.call(rbind, lapply(files_b, read_bed_file))

# join the two tables by matching the name column and discard any rows without a match
tbl <- merge(tbl_a, tbl_b, by = "name", all = FALSE, suffixes = c(".a", ".b"))

# print number of rows in the final table
cat("Number of rows in the final table: ", nrow(tbl), "\n")

# make the plot
r <- ggplot(tbl, aes(x = log10(CRAC_score.a), y = log10(CRAC_score.b))) +
  geom_point(alpha = 0.5, size = 0.5, color = "black") +
  theme_classic(base_size = 60) +
  theme(axis.text = element_text(colour = "black"), axis.title.x = element_blank(),
        axis.title.y = element_blank())

ggsave(commandArgs(trailingOnly = TRUE)[3], plot = r, dpi = 400, width = 10, height = 10)

# print r squared value for the fit between log10(CRAC_score.a) and log10(CRAC_score.b)
fit <- lm(log10(CRAC_score.b) ~ log10(CRAC_score.a), data = tbl)
cat("R squared value for the fit between log10(CRAC_score.a) and log10(CRAC_score.b): ", summary(fit)$r.squared, "\n")
cat("Other fit statistics:\n")
print(summary(fit))

# make a sankey diagram if requested
if (length(commandArgs(trailingOnly = TRUE)) < 4) {
  cat("Sankey diagram not requested.\n")
} else {
  # this code was mostly copied from the ggsankey package documentation
  df <- tbl %>% make_long(VS.a, VS.b)

  s <- ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node,
                        fill = factor(node), label = node)) +
        geom_sankey(flow.alpha = .6, node.color = "gray30", width = 0.3) +
        geom_sankey_label(size = 12, color = "white", fill = "gray40") +
        scale_fill_viridis_d(drop = FALSE) +
        theme_sankey(base_size = 40) +
        labs(x = NULL) +
        theme(legend.position = "none",
              plot.title = element_blank())

  ggsave(commandArgs(trailingOnly = TRUE)[4], plot = s, dpi = 400, width = 10, height = 10)
}

# print entries that were in Q1 of set A but QN of set B and vice versa, where N is the highest quartile
# i.e. the one with the lowest CRAC score.
# these are like outliers in the sankey diagram, there may be something interesting here.

# first, get the max quartile number for set A and set B
N_A <- as.numeric(gsub("Q", "", tail(sort(unique(tbl$VS.a)), 1)))
N_B <- as.numeric(gsub("Q", "", tail(sort(unique(tbl$VS.b)), 1)))

# print rows that were in Q1 of set A but QN of set B
cat("Entries that were in Q1 of set A but Q", N_B, " of set B:\n")
print(tbl[tbl$VS.a == paste0("Q", 1) & tbl$VS.b == paste0("Q", N_B),])

# print rows that were in Q1 of set B but QN of set A
cat("Entries that were in Q1 of set B but Q", N_A, " of set A:\n")
print(tbl[tbl$VS.b == paste0("Q", 1) & tbl$VS.a == paste0("Q", N_A),])