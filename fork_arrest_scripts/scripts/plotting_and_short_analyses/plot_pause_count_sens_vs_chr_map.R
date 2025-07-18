options(bitmapType = "cairo")
require(RIdeogram)

# load data and plot it.
args <- commandArgs(trailingOnly = TRUE)

# load data and set column names
if (length(args) != 4) {
    cat("Usage: Rscript <scriptName>.R fai_file pause_count_vs_coords pause_sens_vs_coords output_svg_file_name\n")
    cat("fai_file: The .fai file associated with the reference genome\n")
    cat(paste0("pause_count_vs_coords: Pause count vs genomic coordinate. A space-separated file with columns ",
               "and column names (comments start with '#'): Chr Start End Value\n"))
    cat(paste0("pause_sens_vs_coords: Expected pause count vs genomic coordinate. A space-separated file with columns ",
               "and column names (comments start with '#'): Chr Start End Value\n"))
    cat("output_svg_file_name: The name of the output SVG file\n")
    cat("NOTE: both files above must contain coordinates on a sufficiently large window size e.g.: 1 kb\n")
    quit()
}

fname_yeast_contig_lens <- args[1]
fname_pause_count_vs_coords <- args[2]
fname_pause_sens_vs_coords <- args[3]
fname_svg <- args[4]

# check that the files above exist
if (!file.exists(fname_yeast_contig_lens)) {
    cat(paste0("File ", fname_yeast_contig_lens, " does not exist"))
    quit()
}
if (!file.exists(fname_pause_count_vs_coords)) {
    cat(paste0("File ", fname_pause_count_vs_coords, " does not exist"))
    quit()
}
if (!file.exists(fname_pause_sens_vs_coords)) {
    cat(paste0("File ", fname_pause_sens_vs_coords, " does not exist"))
    quit()
}

# make the output directory if it does not exist
dir_out <- dirname(fname_svg)
if (!dir.exists(dir_out)) {
    dir.create(dir_out, recursive = TRUE)
}

# load data
table_yeast_contig_lens <- read.table(fname_yeast_contig_lens,
    sep = "\t", header = FALSE, comment.char = "#", stringsAsFactors = F)

table_pause_count_vs_coords <- read.table(fname_pause_count_vs_coords,
    sep = " ", header = FALSE, comment.char = "#", stringsAsFactors = F)

table_pause_sens_vs_coords <- read.table(fname_pause_sens_vs_coords,
    sep = " ", header = FALSE, comment.char = "#", stringsAsFactors = F)

colnames(table_yeast_contig_lens) <- c("Chr", "End", "ignore_1", "ignore_2", "ignore_3")
table_yeast_contig_lens$Start <- 0
colnames(table_pause_count_vs_coords) <- c("Chr", "Start", "End", "Value")
colnames(table_pause_sens_vs_coords) <- c("Chr", "Start", "End", "Value")

table_pause_count_vs_coords$Color <- "fc8d62"

table_yeast_contig_lens <- subset(table_yeast_contig_lens, Chr!= "chrM")
table_yeast_contig_lens <- table_yeast_contig_lens[, c("Chr", "Start", "End") ]

# plot data
ideogram(karyotype = table_yeast_contig_lens,
  overlaid = table_pause_sens_vs_coords,
  label = table_pause_count_vs_coords, label_type = "line",
  output = fname_svg)
