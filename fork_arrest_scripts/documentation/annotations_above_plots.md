# Instructions on how to make annotations above rainplots

NOTE: please read `plotting_workflows.md` before reading this note

## Introduction

You may have seen rainplots with an annotation track above the plot like in the image below (please
do not assume this image was made from real data).

![rainplot](sample_rain_plot_w_annotations.jpg "Rainplot")

In this note, we will discuss how to do this.

## Procedure

First, please prepare a bed file with your annotations.
It is advisable that the bed files are in the BED6 format i.e. tab-separated with no header,
with the columns contig, start, end, name, score, strand with no column header.
Any + feature in the bed file that overlaps a read of interest will appear as a blue line atop the rainplot,
a - feature as a red line, and a '.' feature as a black line ('.' in the strand column means strand not known).

Then, please run
```bash
sbatch annotate_mod_bam_with_bed.sh /path/to/input.bam /path/to/output.bam /path/to/sample_a.bed:sample_a \
  /path/to/sample_b.bed:sample_b
```
to annotate the mod bam file. The input file is a mod BAM file, the output is the same as the input except every
read entry has optional `XT` tags added to them specifying the annotations that overlap the read.
You can specify as many bed files as you want after these two arguments, in the format "sample_a.bed:plot_something",
you have to use the string prefix "plot_" to tell the plotting program that these are annotations that have to be
plotted (replace something with some label that is useful to you).

Then you run `plot_n_reads.sh` or another script as usual with the new output.bam file.
If your annotations are very small (e.g. a few bp long) compared to the total length of the read (e.g. 10 kb),
then you may want to set dpi in `plot_one_read_w_win_or_model_if_needed.R` as a higher value, e.g. 300 or 600
so that the annotations are visible when you zoom in using an image viewer.