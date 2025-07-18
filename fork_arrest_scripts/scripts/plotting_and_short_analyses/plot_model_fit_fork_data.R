#!/usr/bin/env Rscript
# usage Rscript --vanilla scriptName.R inputData outputPicture
require(ggplot2)
require(dplyr)
options(bitmapType = "cairo")

# load data and plot it.

args <- commandArgs(trailingOnly = TRUE)

fit_table <- read.table(args[1], header = FALSE, comment.char = "#")

if(ncol(fit_table) == 11 && length(args) >= 2){

  colnames(fit_table) <- c("dummy", "id", "mean_brdU", "start", "end",
      "win_x", "win_rej", "gof", "lowLim", "upLim", "label")

} else if(ncol(fit_table) == 10 && length(args) >= 2){

  colnames(fit_table) <- c("id", "mean_brdU", "start", "end",
    "win_x", "win_rej", "gof", "lowLim", "upLim", "label")

} else {

  stop(paste("Usage: Rscript ./<script_name>.R plot_data plot.png [true]. ",
             "       plot_data must have 10 or 11 columns",
             "       If the third (optional) argument is 'true', we use publication aesthetics for the plots.",
             ""), sep="\n", call.= FALSE)

}

# set flag for publication aesthetics if the third argument is 'true'
if(length(args) == 3 && args[3] == "true"){
  is_publication_aesthetics <- TRUE
} else {
  is_publication_aesthetics <- FALSE
}

# Remove rows where label is "wrong" or "bad". This is a relic from an older version of the input data.
fit_table <- fit_table[ !(fit_table$label == 'wrong' | fit_table$label == 'bad'), ]

# convert all good_short to good_long in label
fit_table$label <- gsub("good_short", "good_long", fit_table$label)

# An older version of input data had "bad", "good_short", and the newer version doesn't.
# So, the script switches between the two versions suitably in assigning colors and labels.
if("good_short" %in% fit_table$label | "bad" %in% fit_table$label){
  colors <- c("aggData" = "#910c00ff",
      "bad" = "#713c0011",
      "good_long" = "#0000ff11",
      "good_short" = "#0044c011",
      "model" = "#005e00ff")
  labels <- c("All forks average (+- SD)",
                "Bad forks",
                "Good forks >= 30 kb",
                "Good forks < 30 kb",
                "Model")
} else {
  colors <- c("aggData" = "#910c00ff",
      "good_long" = "#0000ff11",
      "model" = "#005e00ff")
  labels <- c("Aggregate forks (mean +- SD)",
               "Individual forks",
               "Model (mean +- SD)")
}

# Drop all rows in fit_table where win_rej is 1
fit_table <- fit_table[fit_table$win_rej != 1, ]

# An older version of input data had accepted and rejected windows and the newer version doesn't.
# So, the script switches between the two versions suitably and plots the data.
if(1 %in% fit_table$win_rej){
  fit_table$win_rej <- as.factor(fit_table$win_rej)
  win_rej_labels <- c("0" = "no", "1" = "yes")

  p <- ggplot(fit_table, aes(x = win_x, y = mean_brdU, group = id)) +
          geom_point(aes(color = label, shape = win_rej), show.legend = TRUE) +
          scale_shape_discrete(name = "Data rejected for model fitting",
              labels = win_rej_labels)
} else if(!is_publication_aesthetics){
  p <- ggplot(fit_table, aes(x = win_x, y = mean_brdU, group = id)) +
          geom_point(aes(color = label), show.legend = TRUE, shape="circle", size=3)
} else {
  p <- ggplot(fit_table, aes(x = win_x, y = mean_brdU, group = id)) +
          geom_point(data = . %>% filter(label != "model"),
                     aes(color = label), show.legend = TRUE, shape="circle", size=3) +
         geom_line(data = . %>% filter(label == "model"),
                   aes(color = label), show.legend = TRUE, linewidth = 2)
}

# add features common to all versions
p <- p + geom_ribbon(aes(ymin = mean_brdU - lowLim, ymax = mean_brdU + upLim,
              fill = label), alpha = .25) +
          xlab("Window coordinate (kb)") +
          ylab("Mean BrdU density (B/B+T)") +
          scale_fill_manual(name = "Legend",
              labels = labels,
              values = colors,
              aesthetics = c("colour", "fill", "line"))

# set publication aesthetics if the flag is set
if(is_publication_aesthetics){
  p <- p + ylab("BrdU density") + xlab("Replisome track coordinate (kb)") +
    theme_classic(base_size = 30) + theme(legend.position="none")
  ggsave(args[2], plot = p, dpi = 400, width = 10, height = 7)
} else{
  # just save the plot
  ggsave(args[2], plot = p, dpi = 400)
}

# print some statistics and plot some statistics
# ----------------------------------------------
# drop rows where label is aggData or model, and group table by win_x, and summarize the number of rows per group
fit_table_stats <- fit_table[!(fit_table$label %in% c("aggData", "model")), ] %>%
  group_by(win_x) %>%
  summarize(n = n())

# print the number of rows per group
print(fit_table_stats, n = nrow(fit_table_stats))

# plot the above data
p <- ggplot(fit_table_stats, aes(x = win_x, y = n)) +
  geom_line(size=2.5) +
  xlab("Window coordinate (kb)") +
  ylab("Number of data points")

# set publication aesthetics if the flag is set
if(is_publication_aesthetics){
  p <- p  + ylim(0, 500) + scale_y_continuous(breaks = c(0,200,400)) +
    theme_classic(base_size = 60) + theme(axis.title.x = element_blank(),
                                          axis.title.y = element_blank(),
                                         panel.background = element_rect(fill='transparent'),
                                         plot.background = element_rect(fill='transparent', color=NA))
  ggsave(paste0(gsub(".png", "", args[2]), "_n_forks.png"), plot = p, dpi = 400, width = 8, height = 5)
} else {
  # just save the plot
  ggsave(paste0(gsub(".png", "", args[2]), "_n_forks.png"), plot = p, dpi = 400)
}

# save the number of forks per window coordinate to a file
write.table(fit_table_stats, paste0(gsub(".png", "", args[2]), "_n_forks.tsv"),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

# print total number of forks in the table
print("Total number of forks: ")
print(fit_table[!(fit_table$label %in% c("aggData", "model")), ] %>%
  group_by(id) %>%
  summarize(n = n()) %>% nrow())

# now, let's calculate sum of residuals
# extract all rows where label is model
model_table <- fit_table[fit_table$label == "model", ]

# drop all rows where label is model or aggData
fit_table <- fit_table[!(fit_table$label %in% c("model", "aggData")), ]

# for each row in fit_table, find the corresponding row in model_table, and subtract mean_brdU in fit_table
# from mean_brdU in model_table
residuals <- fit_table %>%
  mutate(residual = mean_brdU - model_table$mean_brdU[match(win_x, model_table$win_x)])

# print sum of residuals
print("Mean of residuals (comparing each point per fork to model): ")
print(sum(residuals$residual)/length(residuals$residual))

# print sum of squared residuals
print("Mean of square of residuals (comparing each point per fork to model): ")
print(sum(residuals$residual^2)/length(residuals$residual))