require(ggplot2)
require(dplyr)
options(bitmapType = "cairo")

# goal of program
# ===============
# plot bar plot of expected and observed pause counts versus some division of a genomic feature
# (e.g.: quartiles of tRNA grouped according to transcription rates)

# typical usage
# ==============
# < some_file.txt Rscript <script_name>.R /path/to/output_plot.png xlabel ylabel width height rel_text_size_x \
#     rel_text_size_y rel_title_size_x rel_title_size_y width_whisker_rel_bar remove_columns division_conversion_rule \
#     ylabel_ratio ylims_ratio sig_threshold
# NOTE: \ means that the command continues on the next line. There should be no space after the \.
# some_file.txt: a file with columns separated by one or more space or tabs,
#      with column names expected, sd, observed, division, coverage, coverage_sd, expected_sens, expected_sens_sd.
#      Only observed and division are required. Any other graph is produced if other columns are suitably present.
#      Column meanings are as follows:
#      expected is the expected pause count assuming uniform pause probability per base across each fork,
#      sd is the standard deviation of the expected pause count,
#      observed is the observed pause count,
#      division is a label for the division of the genomic feature
#        (e.g.: quartiles of tRNA grouped according to transcription rates),
#      coverage is the fraction of the genome covered by the ROI,
#      coverage_sd is the standard deviation of the coverage,
#      coverage_T is the fraction of genomic thymidines covered by the ROI (the column could have been restricted to one
#        or both strands, we leave it to the user/the calling script to provide appropriate input here),
#      coverage_T_sd is the standard deviation of the coverage_T,
#      expected_sens is the expected pause count from sensitivity calculations, and
#      expected_sens_sd is the standard deviation of the expected pause count from sensitivity calculations.
# output_plot_name.png: output plot name, must end in .png.
#                       The program will also create a plot of the ratio of observed to expected sensitivity if
#                       expected_sens and expected_sens_sd are present in the input file.
# Henceforth, every argument is optional.
# xlabel: label of x axis, defaults to "Rank"
# ylabel: label of y axis, defaults to "Fraction of global pause count (%)".
#         NOTE: ylabel of the ratio plot is set separately later in the argument list.
# width: width of plot in inches, defaults to 5
# height: height of plot in inches, defaults to 5
# rel_text_size_x: relative size of x axis text, defaults to 1
# rel_text_size_y: relative size of y axis text, defaults to 1
# rel_title_size_x: relative size of x axis title, defaults to 1
# rel_title_size_y: relative size of y axis title, defaults to 1
# width_whisker_rel_bar: width of the error bar whiskers relative to the width of the bars, defaults to 0.2
# remove_columns: a comma separated list of columns to remove from the input file before plotting, defaults to none
# division_conversion_rule: rules to convert divisions to numeric, default to no conversion,
#                           leave it blank "" if not needed.
#                           Format is a string with three entries separated by commas: zero_bin,is_increasing,bin_size.
#                           See the details below.
# ylabel_ratio: label of y axis for the ratio plot, defaults to "Observed:expected count".
# ylims_ratio: y limits for the ratio plot in the format num1,num2; defaults to 0,2.
# sig_threshold: threshold for significance, defaults to 3 i.e. > or < 3 s.d.
#                - You can specify two values separated by a comma here, then the higher value and lower values
#                  will be marked with a '*' and a '†' respectively.
#                  You have to specify the higher value first and the lower value second.
#                - You can specify three items here, then the first two values will be treated as above, and the
#                  only allowed value for the third entry is "no_depletion", in which case any significant result
#                  due to depletion i.e. due to the observed value being lower than the expected value will be ignored.
#                  This is useful if we are only interested in calculating enrichment and depletion is not of interest.

# Details of division conversion rule
# ===================================
# Our division column looks like "LS_12", "LS_13" etc. or like "VS_1_LS_13" etc. where VS stands for value split and
# LS stands for length split. So LS_13 means 13 units of length, LS_14 means 14 units of length etc.
# So, if we know the zero position, the bin size, and whether the division is increasing or decreasing, we can convert
# the division to numeric. For example, if the zero position is 12, the bin size is 1 bp,
# and the division is increasing, then LS_12 will be converted to 0, LS_13 will be converted to 1 bp,
# LS_14 will be converted to 2 etc.
# So we need to know three pieces of information: zero_bin, is_increasing, bin_size.
# For the example above, zero_bin is 12, is_increasing is TRUE, and bin_size is 1.
# More examples
# -------------
# Bins: VS_1_LS_1, VS_1_LS_2, VS_1_LS_3, VS_1_LS_4, VS_1_LS_5, VS_1_LS_6, VS_1_LS_7, VS_1_LS_8, VS_1_LS_9, VS_1_LS_10
# zero_bin: 5, is_increasing: TRUE, bin_size: 100
# Converted to: -400, -300, -200, -100, 0, 100, 200, 300, 400, 500
# Bins: VS_1_LS_1, VS_1_LS_2, VS_1_LS_3, VS_1_LS_4, VS_1_LS_5, VS_1_LS_6, VS_1_LS_7, VS_1_LS_8, VS_1_LS_9, VS_1_LS_10
# zero_bin: 5, is_increasing: FALSE, bin_size: 1000
# Converted to: 4000, 3000, 2000, 1000, 0, -1000, -2000, -3000, -4000, -5000
# Bins: VS_1_LS_1, VS_1_LS_2, VS_2_LS_3
# zero_bin: 5, is_increasing: FALSE, bin_size: 1000
# This would throw an error because the division prefixes are not the same i.e. VS_1 and VS_2 are different.

# load command line arguments
args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 1){
    print("Error: must supply at least one argument. The other arguments are optional.")
    stop(cat("Usage: < some_file.txt Rscript <script_name>.R /path/to/output_plot.png xlabel ylabel width height
                       rel_text_size_x rel_text_size_y rel_title_size_x rel_title_size_y width_whisker_rel_bar
                       remove_columns division_conversion_rule ylabel_ratio ylims_ratio sig_threshold"))
}

# set default values for input arguments
x_label <- "Rank"
y_label <- "Fraction of global pause count (%)"
width <- 5
height <- 5
rel_text_size_x <- 1
rel_text_size_y <- 1
rel_title_size_x <- 1
rel_title_size_y <- 1
width_whisker_rel_bar <- 0.2
ylabel_ratio <- "Observed:expected count"
y_lim_lower_ratio <- 0
y_lim_upper_ratio <- 2

# set to user supplied labels if available
if(length(args) >= 2 && args[2] != ""){
    x_label <- args[2]
}

if(length(args) >= 3 && args[3] != ""){
    y_label <- args[3]
}

# set width and height if user supplied
if(length(args) >= 4 && args[4] != ""){
    width <- as.numeric(args[4])
}

if(length(args) >= 5 && args[5] != ""){
    height <- as.numeric(args[5])
}

# set relative sizing variables if supplied
if(length(args) >= 6 && args[6] != ""){
    rel_text_size_x <- as.numeric(args[6])
}
if(length(args) >= 7 && args[7] != ""){
    rel_text_size_y <- as.numeric(args[7])
}
if(length(args) >= 8 && args[8] != ""){
    rel_title_size_x <- as.numeric(args[8])
}
if(length(args) >= 9 && args[9] != ""){
    rel_title_size_y <- as.numeric(args[9])
}
if(length(args) >= 10 && args[10] != ""){
    width_whisker_rel_bar <- as.numeric(args[10])
}
if(length(args) >= 11 && args[11] != ""){
    remove_columns <- unlist(strsplit(args[11], ","))
}
is_division_conversion_rule <- FALSE
if(length(args) >= 12 && args[12] != ""){
    division_conversion_rule <- unlist(strsplit(args[12], ","))
    if(length(division_conversion_rule) != 3){
      is_division_conversion_rule <- FALSE
    } else {

      is_division_conversion_rule <- TRUE

      zero_bin <- as.numeric(division_conversion_rule[1])
      is_increasing <- as.logical(division_conversion_rule[2])
      bin_size <- as.numeric(division_conversion_rule[3])

      # raise an error if any of the values are not numeric
      if(any(is.na(c(zero_bin, is_increasing, bin_size)))){
          stop(cat("Error: division_conversion_rule must have numeric values."))
      }

      # raise an error if bin_size is not positive
      if(bin_size <= 0){
          stop(cat("Error: bin_size must be positive."))
      }

    }
}
if(length(args) >= 13 && args[13] != ""){
    ylabel_ratio <- args[13]
}
if(length(args) >= 14 && args[14] != ""){
    ylims_ratio <- unlist(strsplit(args[14], ","))
    if(length(ylims_ratio) != 2){
        stop(cat("Error: ylims_ratio must have two values."))
    }
    y_lim_lower_ratio <- as.numeric(ylims_ratio[1])
    y_lim_upper_ratio <- as.numeric(ylims_ratio[2])
}

if(length(args) >= 15 && args[15] != ""){

  sig_list_args <- unlist(strsplit(args[15], ","))

  # first check if the user has supplied many values separated by a comma
  if (length(sig_list_args) == 2 || length(sig_list_args) == 3){

    sig_threshold <- as.numeric(sig_list_args[1])
    sig_threshold_lower <- as.numeric(sig_list_args[2])

    # complain if the user has supplied a lower threshold that is higher than the higher threshold
    if (sig_threshold <= sig_threshold_lower){
      stop(cat("Error: the lower threshold must be less than the higher threshold."))
    }

    # check if the user has supplied a third argument and if it is "no_depletion" set the flag
    no_calculate_depletion_sig <- (length(sig_list_args) == 3 && sig_list_args[3] == "no_depletion")

    # complain if the third argument is not "no_depletion"
    if (length(sig_list_args) == 3 && !no_calculate_depletion_sig){
      stop(cat("Error: the third argument must be 'no_depletion' if supplied."))
    }

  } else {
    sig_threshold <- as.numeric(sig_list_args[1])
    sig_threshold_lower <- Inf
    no_calculate_depletion_sig <- FALSE
  }
} else{
    sig_threshold <- 3
    sig_threshold_lower <- Inf
    no_calculate_depletion_sig <- FALSE
}

# input values in the array below
try_reading_file <- try(read.table(file("stdin"), comment.char = "#", header = TRUE))
if(class(try_reading_file) == "try-error"){
    quit()
} else {
    df2 <- try_reading_file
}

# remove columns if requested
if(exists("remove_columns")){
    df2 <- df2[, !colnames(df2) %in% remove_columns]
}

# check that the columns division, expected, sd, observed, and coverage are present
if(!all(c("division", "observed") %in% colnames(df2))){
    stop(cat("Error: input file must have columns division, observed at the minimum."))
}

# check whether the division column is not numeric
is_non_numeric_division <- is.character(df2$division)

# convert division if requested
if(is_division_conversion_rule){

  if(is_non_numeric_division){

    # set the non-numeric division flag to FALSE
    is_non_numeric_division <- FALSE

    # split division by _ and retain just the last element
    df2$division_num <- as.numeric(sapply(strsplit(df2$division, "_"), function(x) as.numeric(x[length(x)])))

    # get the division prefix
    df2$division_prefix <- sapply(strsplit(df2$division, "_"), function(x) paste(x[-length(x)], collapse="_"))

    # ensure that all the division prefixes are the same
    if(length(unique(df2$division_prefix)) > 1){
      stop(cat("Error: division prefixes are not the same."))
    }

    # convert division into an actual length scale
    if(is_increasing){
      df2$division <- (df2$division_num - zero_bin) * bin_size
    } else{
      df2$division <- (zero_bin - df2$division_num) * bin_size
    }

  } else {
    stop(cat("Error: division_conversion_rule is only applicable to non-numeric divisions."))
  }
}

# if division is a string, ensure that the input order of the divisions is preserved
if(is_non_numeric_division){
  division_list <- unique(df2$division)
  df2 <- within(df2, division <- factor(division, levels=division_list))
} else{
  xlim <- c(min(df2$division), max(df2$division))
  bin_size <- abs(diff(xlim))/length(df2$division)
  xlim <- c(xlim[1] - bin_size/2, xlim[2] + bin_size/2)
}


# add a column to indicate whether the division is significant
if(all(c("expected_sens", "expected_sens_sd") %in% colnames(df2))){

  if (no_calculate_depletion_sig){
    df2$sig_level <- (df2$observed - df2$expected_sens)/df2$expected_sens_sd
  } else {
    df2$sig_level <- abs(df2$observed - df2$expected_sens)/df2$expected_sens_sd
  }

  df2 <- df2 %>% mutate(significant = case_when(
      sig_level >= sig_threshold ~ "yes",
      sig_level >= sig_threshold_lower  ~ "yes_but_lower",
      TRUE ~ "no"
  ))
}

# scale the different fractions by a factor 1/100.
# This is because the fractions are stored in basis points (1 basis point = 0.01%),
# and we want to plot them in percentage points.
for (col in c("observed", "expected", "expected_sens", "coverage", "expected_sens_sd", "sd", "coverage_sd",
              "coverage_T", "coverage_T_sd")){
    if(col %in% colnames(df2)){
        df2[,col] <- df2[,col] / 100
    }
}

# set y limits
y_lim <- max(df2$observed) * 1.1
is_sensitivity <- FALSE
if(all(c("expected", "sd") %in% colnames(df2))){
    y_lim <- max(y_lim, max(df2$expected + df2$sd) * 1.1)
}
if(all(c("coverage", "coverage_sd") %in% colnames(df2))){
    y_lim <- max(y_lim, max(df2$coverage + df2$coverage_sd) * 1.1)
}
if(all(c("coverage_T", "coverage_T_sd") %in% colnames(df2))){
    y_lim <- max(y_lim, max(df2$coverage_T + df2$coverage_T_sd) * 1.1)
}
if(all(c("expected_sens", "expected_sens_sd") %in% colnames(df2))){
    is_sensitivity <- TRUE
    y_lim <- max(y_lim, max(df2$expected_sens + df2$expected_sens_sd) * 1.1)
}

# Bar plots with optional tracks depending on input columns
p <- ggplot(df2, aes(x=division, y=observed)) +
     geom_bar(stat="identity", color="#DDCC77", fill="#DDCC77",
              position=position_dodge()) +
     if(all(c("expected", "sd") %in% colnames(df2))){
         list(geom_errorbar(aes(ymin=expected-sd, ymax=expected+sd), width=width_whisker_rel_bar,
                       position=position_dodge(.9), size=1, color="#000000"),
              geom_point(aes(y=expected), position=position_dodge(.9), size=1, color="#000000"))
     }

p <- p +
     if(all(c("coverage", "coverage_sd") %in% colnames(df2))){
         list(geom_errorbar(aes(ymin=coverage-coverage_sd, ymax=coverage+coverage_sd), width=width_whisker_rel_bar,
                       position=position_dodge(.9), size=1, color="#332288", linetype="dotted"),
              geom_point(aes(y=coverage), position=position_dodge(.9), size=1, color="#332288"))
     }

p <- p +
     if(all(c("coverage_T", "coverage_T_sd") %in% colnames(df2))){
         list(geom_errorbar(aes(ymin=coverage_T-coverage_T_sd, ymax=coverage_T+coverage_T_sd),
                            width=width_whisker_rel_bar, position=position_dodge(.9), size=1, color="#332288"),
              geom_point(aes(y=coverage_T), position=position_dodge(.9), size=1, color="#332288", shape=0))
     }

p <- p +
     if(all(c("expected_sens", "expected_sens_sd") %in% colnames(df2))){
         list(geom_errorbar(aes(ymin=expected_sens-expected_sens_sd, ymax=expected_sens+expected_sens_sd),
                       width=width_whisker_rel_bar, position=position_dodge(.9), size=1, color="#117733"),
              geom_point(aes(y=expected_sens), position=position_dodge(.9), size=1, color="#117733"),
              geom_text(aes(label = case_when(significant == "yes" ~ "*",
                                        significant == "yes_but_lower" ~ "†",
                                        significant == "no" ~ "")),
                       position = position_dodge(width = .9), vjust = -0.2, size = 20 / .pt))
     }

# Finished bar plot
p <- p + labs(x=x_label, y=y_label)

if(is_non_numeric_division){
  p <- p + coord_cartesian(ylim=c(0, y_lim))
} else{
  p <- p + coord_cartesian(xlim=xlim, ylim=c(0, y_lim))
}

p <- p + theme_classic() +
         scale_fill_manual(values= '#999999')

if(is_non_numeric_division){
  p <- p + scale_x_discrete(guide = guide_axis(angle = 45))
}

p <- p + theme_classic(base_size = 22) +
         theme(axis.text.x=element_text(size=rel(rel_text_size_x)),
               axis.text.y=element_text(size=rel(rel_text_size_y)),
               axis.title.x=element_text(size=rel(rel_title_size_x)),
               axis.title.y=element_text(size=rel(rel_title_size_y)))

ggsave(args[1], plot = p, dpi = 400, width = width, height = height)

# if sensitivity is present, then plot a ratio of observed to expected sensitivity
if(is_sensitivity){
    p <- ggplot(df2, aes(x=division, y=observed/expected_sens)) +
         geom_bar(stat="identity", color="#DDCC77", fill="#DDCC77",
                  position=position_dodge()) +
         geom_text(aes(label = case_when(significant == "yes" ~ "*",
                                        significant == "yes_but_lower" ~ "†",
                                        significant == "no" ~ "")),
                       position = position_dodge(width = .9), vjust = -0.2, size = 20 / .pt)

    # Finished bar plot
    p <- p + labs(x=x_label, y=ylabel_ratio)

    if(is_non_numeric_division){
      p <- p + coord_cartesian(ylim=c(y_lim_lower_ratio, y_lim_upper_ratio))
    } else{
      p <- p + coord_cartesian(xlim=xlim, ylim=c(y_lim_lower_ratio, y_lim_upper_ratio))
    }

    p <- p + theme_classic() +
             scale_fill_manual(values= '#999999')

    if(is_non_numeric_division){
      p <- p + scale_x_discrete(guide = guide_axis(angle = 45))
    }

    p <- p + theme_classic(base_size = 22) +
             theme(axis.text.x=element_text(size=rel(rel_text_size_x)),
                   axis.text.y=element_text(size=rel(rel_text_size_y)),
                   axis.title.x=element_text(size=rel(rel_title_size_x)),
                   axis.title.y=element_text(size=rel(rel_title_size_y)))

    ggsave(sub(".png", "_ratio.png", args[1]), plot = p, dpi = 400,
                  width = width, height = height)
}