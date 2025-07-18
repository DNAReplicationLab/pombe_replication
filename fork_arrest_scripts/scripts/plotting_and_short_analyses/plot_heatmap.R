#!/usr/bin/env Rscript
# Written by Isabel Díez Santos and Sathish Thiyagarajan

require(ggplot2)
require(tidyverse)
require(viridis)
options(bitmapType = "cairo")

# Function to plot heatmap
plot_heatmap <- function(input_file, output_file, highlight_bed_file = NULL, rDNA_or_not = FALSE,
                         align_by = FALSE, x_limits = NULL, y_limits = NULL, flip_by_strand = FALSE) {

  # check that flip_by_strand is only set to TRUE if align_by is also set to TRUE
  if (flip_by_strand && !align_by) {
      stop("flip_by_strand can only be set to TRUE if align_by is also set to TRUE")
  }

  # Read the data from the input file
  data <- read.table(input_file, comment = "#", header=TRUE)

  # Check if data contains expected columns
  required_cols <- c("read_id", "start", "end", "mod_qual")
  if (!all(required_cols %in% colnames(data))) {
    stop("Input file must contain the following columns: read_id, start, end, mod_qual and optionally alpha")
  }

  # assign column types
  data$read_id <- as.character(data$read_id)
  data$start <- as.double(data$start)
  data$end <- as.double(data$end)
  data$mod_qual <- as.double(data$mod_qual)
  if ("alpha" %in% colnames(data)) {
      data$alpha <- as.double(data$alpha)
  }

  # if a highlight_bed_file is provided and align by is requested,
  # read the file (retain only the first row for each read id i.e. we can only afford one highlight per read id),
  # drop all rows where read id is not in the highlight bed file,
  # and subtract the highlight start from start and end to perform the alignment
  if (!is.null(highlight_bed_file) && align_by) {

    highlight_data <- read.table(highlight_bed_file, header = FALSE, comment.char = "#")
    colnames(highlight_data)[1:6] <- c("contig", "start", "end", "read_id", "score", "strand")
    highlight_data$start <- as.double(highlight_data$start)

    highlight_data <- highlight_data %>% group_by(read_id) %>% slice(1) %>% ungroup()

    data <- data %>% filter(read_id %in% highlight_data$read_id)

    if(nrow(highlight_data) == 0) {
      stop("No highlights found in the input file")
    } else {
      for (i in 1:nrow(highlight_data)) {
        rows <- data$read_id == highlight_data$read_id[i]
        data$start[rows] <- data$start[rows] - highlight_data$start[i]
        data$end[rows] <- data$end[rows] - highlight_data$start[i]
        if (flip_by_strand && highlight_data$strand[i] == "-") {
          data$start[rows] <- -data$start[rows]
          data$end[rows] <- -data$end[rows]
        }
      }
    }
  }

  # get the mean value of mod_qual per read_id
  data_max <- data %>% group_by(read_id) %>%
      summarise(mean_mod_qual = mean(mod_qual, na.rm = TRUE)) %>%
      ungroup()

  # retain those reads where the mean value of mod_qual is greater than 0.05
  data_max <- data_max %>% filter(mean_mod_qual > 0.05)
  data <- data %>% filter(read_id %in% data_max$read_id)

  # Remove rows with NA values in critical columns
  data <- data %>% drop_na(read_id, start, mod_qual)

  # calculate width and center of each tile
  data$width <- (data$end - data$start)
  data$center <- (data$start + data$end)/2

  # Determine the first start position and last end position for each read_id
  first_start_last_end <- data %>%
    group_by(read_id) %>%
    summarise(first_start = min(start, na.rm = TRUE),
              last_end = max(end, na.rm = TRUE)) %>%
    arrange(first_start)

  # set up a list of y positions for each unique read_id
  unique_read_ids <- first_start_last_end$read_id
  read_id_index <- seq(from = 1,to = length(unique_read_ids), by = 1)
  last_end_positions_by_read_id <- first_start_last_end$last_end
  if(length(unique_read_ids) == 0) {
    stop("No reads with mod_qual > 0.1 found in the input file")
  } else {
    for (read_i in 1:length(unique_read_ids)) {
      row_j <- 1
      while(!is.na(last_end_positions_by_read_id[row_j]) && row_j < read_i &&
                first_start_last_end$first_start[read_i] <= last_end_positions_by_read_id[row_j] + 5000) {
        row_j <- row_j + 1;
      }
      if (row_j < read_i) {
        read_id_index[read_i] <- row_j;
        last_end_positions_by_read_id[row_j] <- last_end_positions_by_read_id[read_i];
        last_end_positions_by_read_id[read_i] <- NA;
      }
    }
  }

  # add a read_id_index column to the data using unique_read_ids and read_id_index
  data <- data %>% mutate(read_id_index = read_id_index[match(read_id, unique_read_ids)])

  # one line function to round a number down to the nearest multiple of 20000
  round_to_20000 <- function(x) 20000 * floor(x / 20000)

  # Determine the explicit range for the x-axis based on the actual data or the user-provided limits.
  # Throw out data if necessary.
  x_min <- min(data$start, na.rm = TRUE)
  x_max <- max(data$end, na.rm = TRUE)
  if (!is.null(x_limits)) {
    x_min <- x_limits[1]
    x_max <- x_limits[2]
    breaks <- x_limits[3]
    x_breaks <- seq(x_min, x_max, breaks)
    data <- data %>% filter(center >= x_min & center <= x_max)
  } else {
    x_breaks <- seq(round_to_20000(x_min), round_to_20000(x_max), 20000)
  }
  x_break_labels <- x_breaks / 1000  # Convert to kb

  # if y_limits are provided, filter the data
  if (!is.null(y_limits)) {
      data <- data %>% filter(read_id_index >= y_limits[1] & read_id_index <= y_limits[2])
  }

  # if a highlight_bed_file is provided and align by is not requested,
  # read the file and add a highlight column to the data
  if (!is.null(highlight_bed_file) && !align_by) {
    highlight_data <- read.table(highlight_bed_file, header = FALSE, comment.char = "#")
    colnames(highlight_data)[1:6] <- c("contig", "start", "end", "read_id", "score", "strand")
    highlight_data$start <- as.double(highlight_data$start)

    data$highlight <- FALSE
    data$strand <- "."
    if(nrow(highlight_data) == 0) {
      stop("No highlights found in the input file")
    } else {
      for (i in 1:nrow(highlight_data)) {
        # locate rows in data that overlap with the highlight region
        rows <- data$start <= highlight_data$start[i] & data$end > highlight_data$start[i] &
          data$read_id == highlight_data$read_id[i]
        data$highlight[rows] <- TRUE
        data$strand[rows] <- highlight_data$strand[i]
      }
    }
  }

  # Plot the heatmap with customizations
  p <- ggplot(data, aes(x = center, y = read_id_index, width = width, height = 1, fill = mod_qual))

  if ("alpha" %in% colnames(data)) {
    p <- p + geom_tile(aes(alpha = alpha))
  } else {
    p <- p + geom_tile()
  }

  if (!is.null(highlight_bed_file) && !align_by) {
    p <- p +  geom_point(data = data[data$highlight & data$strand == "+",], aes(x = center, y = read_id_index + 1),
                 shape = 25, size = 2, fill = "red") +
              geom_point(data = data[data$highlight & data$strand == "-",], aes(x = center, y = read_id_index - 1),
                 shape = 24, size = 2, fill = "red")
  }

  if (rDNA_or_not && !align_by) {
    for (i in seq(from=8855, to=200732, by=9137)){
      p <- p + geom_vline(xintercept = i, linetype="solid", color = "grey", size=1, alpha=0.3);
    }
  }

  p <- p + scale_fill_viridis(name = "BrdU fraction", option = "D", breaks = c(0, 0.25, 0.5, 0.75, 1),
                       limits = c(0, 1), na.value="white") +  # Explicitly set breaks and limits
    labs(x = "Reference coordinate (kb)", y = "Read index") +
    scale_x_continuous(breaks = x_breaks, limits = c(x_min, x_max), expand = c(0, 0), labels = x_break_labels) +
    scale_y_discrete(expand = c(0, 0)) +  # Remove extra space on y-axis
    theme_bw() +  # Use theme_bw for a clean background
    theme(panel.grid = element_blank(), axis.text.x = element_text(size = 14),
          legend.title = element_text(size = 14), legend.text = element_text(size = 12))  # Remove gridlines

  # Save the plot to a file
  ggsave(paste0(output_file, "_heatmap.png"), p, width = 10, height = 8, device ="png", dpi = 300, create.dir = TRUE)

  # make a list of read_ids and the corresponding read_id_index
  read_ids_and_index <- list(read_ids = unique_read_ids, read_id_index = read_id_index,
                             start = first_start_last_end$first_start, end = first_start_last_end$last_end)

  # Write the list of unique read_ids and corresponding indices to a file
  write.table(read_ids_and_index, file = paste0(output_file_prefix, "_read_ids.txt"),
              row.names = FALSE, col.names = TRUE, quote=FALSE, sep = "\t")

  # Write the data to a file
  write.table(data, file = paste0(output_file, "_heatmap_data.txt"),
              row.names = FALSE, col.names = TRUE, quote=FALSE, sep = "\t")

  # Make an average across the Y axis and write to a file
  # - first, only retain data between the first start and last end positions
  # - calculate a mean interval size
  # - round the center to the nearest multiple of the mean interval size
  # - average the mod_qual values for each center_rounded
  # - write to a file retaining only the center_rounded, n, and avg_mod_qual columns
  data <- data %>% filter(center >= x_min & center <= x_max)
  mean_interval_size <- mean(data$end - data$start)
  data$center_rounded <- mean_interval_size * floor(data$center / mean_interval_size)
  data_avg <- data %>% group_by(center_rounded) %>%
      summarise(avg_mod_qual = mean(mod_qual, na.rm = TRUE), n = n()) %>%
      ungroup()
  write.table(data_avg %>% select(center_rounded, n, avg_mod_qual), file = paste0(output_file, "_heatmap_data_avg.txt"),
              row.names = FALSE, col.names = TRUE, quote=FALSE, sep = "\t")

}

# Main script
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  print("Usage: Rscript plot_heatmap.R <input_file> <output_file_prefix> <highlight_bed_file_optional> \ ")
  print("       <rDNA_optional> <align_by_optional> <x_limits_optional> <y_limits_optional>")
  print(" \ means the command continues on the next line")
  print(" optional means the argument is optional, and you can skip it by providing an empty string.")
  print(" if you want to specify an optional argument, you must specify all the previous optional arguments.")
  print(" you can skip all optional arguments by providing empty strings in all places or just not providing them.")
  print("IMPORTANT NOTE: we ignore contigs throughout this script, so if you provide files with reads on multiple")
  print("                contigs, the X values per read must be interpreted as on the read\'s reference contig. ")
  print("input_file: path to the input file contains read_id, start, end, mod_qual columns, tab-separated with header.")
  print("            you can set this to - to read from stdin. Can optionally contain a column called alpha, ")
  print("            which is used to set the transparency of the tiles (0-1) where 0 is transparent and 1 is opaque.")
  print("output_file_prefix: prefix for the output files e.g. if you supply /a/b/c/output then a few files are written")
  print("                    with their names starting with /a/b/c/output")
  print("highlight_bed_file_optional: path to the highlight bed file, if provided, the regions in the bed file if they")
  print("                             fall within with the regions in the input file, they are highlighted in the")
  print("                             plot. The bed file should have at least 6 columns: contig, start, end, read_id, ")
  print("                             (unused) score, strand. + and - strands are plotted differently. default unused.")
  print("                             You can use a pauseFile here after converting it to bed using the ")
  print("                             convert_pause_file_to_bed.py script (with appropriate options to throw out ")
  print("                             filtered pauses and to convert L/R forks to -/+).")
  print("                             NOTE: we only use the start column from the bed file because we do not want to ")
  print("                                   highlight multiple windows along the read.")
  print("rDNA_optional: if set to rDNA, guidelines are plotted for the rRFB according to our coordinate system, ")
  print("               default unused. If you don't want to use this, do not supply it or set it to anything else.")
  print("align_by_optional: set to 'yes' or 'yes_and_flip_by_strand'; default 'no'.")
  print("                   if set to 'yes', then only reads in the highlight_bed_file are retained in the heatmap,")
  print("                   and the positions are aligned by the start of the highlight region.")
  print("                   if set to 'yes_and_flip_by_strand', then the reads are flipped by strand in the ")
  print("                   highlight_bed_file after alignment. If you\'ve converted pauseFile to bed and set ")
  print("                   the strand column by converting L/R forks to -/+ then you should use this option if you ")
  print("                   want to align by fork direction.")
  print("                   NOTE: if multiple lines correspond to the same read_id in the input file, then only the ")
  print("                         one that appears first in the bed file is retained.")
  print("                   NOTE: depending on how the window boundaries overlap with the highlight bed start ")
  print("                         positions, you may get jitter in the plot when you try to align.")
  print("x_limits_optional: if set to num1,num2,num3 then the x-axis limits are set to num1 and num2 and the breaks")
  print("                   are set to num3. default unused. i.e. the x-axis limits etc. are determined automatically.")
  print("y_limits_optional: if set to num1,num2 then only reads with row index between num1 and num2 are plotted.")
  print("                   default unused. i.e. all reads are plotted. NOTE: we do not know which read is assigned ")
  print("                   which index, so you may need to experiment with this option to get the desired reads.")
  print("                   We recommend plotting the heatmap without this option first and see how it looks, and then")
  print("                   use this option (and perhaps the x_limits_optional) to zoom in on a region of interest. ")
}

input_file <- args[1]
output_file_prefix <- args[2]

# set the optional arguments
if (length(args) >= 3 && args[3] != "" && args[3] != "/dev/null") {
  highlight_bed_file <- args[3]
} else {
  highlight_bed_file <- NULL
}
if (length(args) >= 4 && args[4] == "rDNA") {
  rDNA_or_not <- TRUE
} else {
  rDNA_or_not <- FALSE
}
if (length(args) >= 5) {
  if (args[5] == "yes"){
    align_by <- TRUE
    flip_by_strand <- FALSE
  } else if (args[5] == "yes_and_flip_by_strand") {
    align_by <- TRUE
    flip_by_strand <- TRUE
  } else {
    align_by <- FALSE
    flip_by_strand <- FALSE
  }
} else {
  align_by <- FALSE
  flip_by_strand <- FALSE
}
if (length(args) >= 6 && args[6] != "") {
  x_limits <- as.numeric(strsplit(args[6], ",")[[1]])
} else {
  x_limits <- NULL
}
if (length(args) >= 7 && args[7] != "") {
  y_limits <- as.numeric(strsplit(args[7], ",")[[1]])
} else {
  y_limits <- NULL
}

# if input file is - then read from stdin
if (input_file == "-") {
  input_file <- file("stdin")
}

plot_heatmap(input_file, output_file_prefix, highlight_bed_file, rDNA_or_not, align_by, x_limits, y_limits,
             flip_by_strand)