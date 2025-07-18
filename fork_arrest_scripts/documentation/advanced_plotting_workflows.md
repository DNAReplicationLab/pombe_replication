Please read `plotting_workflows.md` before reading this note.
This note goes into more details on how rain plots are made; please use
the commands here only if the commands in the earlier note do not do what you want.

## Basic plotting workflow

We describe the most basic and direct way to make rain plots.
Do not use this workflow unless you require fine control over the plotting process.

### Main plotting script of analogue modification data

The basic script to make rain plots is `plot_one_read_w_win_or_model_if_needed.R`.
The script is never used directly, but is called by some other script.

Usage:

```shell
Rscript ./plot_one_read_w_win_or_model_if_needed.R $plot_data $plot_filename $plot_elements_this_is_optional
```

Output plot is `$plot_filename` in png format.
The optional plot elements input will be explained later on.

Data for the plot is in the file `$plot_data`.
It is space or tab-separated with five columns but no headers.
The columns are id, start, end, val, and label.
`start` and `end` are coordinates on the reference genome.
`val` is a number between 0 and 1.
`label` can be `rawDetect`, `winDetect`, or `model`.

So, an input file looks like the following (this is fake data).

```text
id1 1 2 0.1 rawDetect
id1 2 3 0.2 rawDetect
id1 3 4 0.3 rawDetect
id1 4 5 0.4 rawDetect
id1 1 3 0.2 winDetect
id1 3 5 0.4 winDetect
id1 1 2 0.3 model
id1 2 3 0.3 model
id1 3 4 0.3 model
id1 4 5 0.3 model
```

The interpretation is that analogue modification probability on the molecule
with read id 'id1' at positions 1, ..., 4 are 0.1, ..., 0.4.
Windowing has resulted in values of 0.2, 0.4 in the windows 1-3 and 3-5 respectively.
A model has predicted values of 0.3 everywhere.
These three types of data are plotted with three different representations in the figure.

NOTE: all of these are wrong values. They are just used for illustration.
Do not be bothered that the windowing is wrong or that the model is poor.
The plotting script does not check these values to see if they make sense.

One can introduce annotations in the plot using the file `$plot_elements_this_is_optional`.
As the name implies, this step is optional; just omit the command line argument to avoid this.
The file has three necessary columns and one optional column separated by spaces and no header.
The necessary columns are start, end, label; the optional column is line style.
start and end are positions on the reference genome.
start must always be less than or equal to end i.e. don't switch them for left/right fork.
Label can be 'origin', 'leftFork', 'rightFork', 'termination', or 'pause'.
Line style is a (binary) column for line style in plot for a fork feature.
Set it to zero or one if you want to distinguish between fork features.

As stated earlier, the script is rarely run by itself, because our working style is that
we run other commands to construct a plot data and/or a plot annotation file on the fly,
and then call this plotting script. This is because our data is stored in other formats
and has to first be recast into this plot data per read id format before this script
can be run.

### Preparing input data for plotting and feeding it to the plotting script 

Mod bam file, read id, contig, start, end, model parameters, and fasta file of the reference genome are
the inputs needed for the snippet below.

```shell
{ 
	# get raw data
	python get_raw_data_from_modBAM.py $mod_bam $read_id chrXVI 917920 948066 |\
	  awk '{print $1 " " $2 " " $2+1 " " $3 " rawDetect"}';
	  
	# window it
	python get_raw_data_from_modBAM.py $mod_bam $read_id chrXVI 917920 948066 |\
	  sed '1i\detectIndex\tposOnRef\tval' |\
	  python get_mean_brdU_window.py --window 300 --thres 0.5 |\
	  awk '{print $1 " " $3 " " $4 " " $2 " winDetect"}';
	  
	# add model curve
	python get_model_values_at_raw_data_coords.py sigmoidal chrXVI 917920 948066\
          0.7_0.1_940000_3000 $fasta_file | awk '{print $1 " " $1+1 " " $2 " model"}' |\
          sed "s/^/$read_id /";
        
} > plotting_and_short_analyses/plot_data

# plot data
cd plotting_and_short_analyses;
Rscript ./plot_one_read_w_win_or_model_if_needed.R plot_data plot.png
```

If only read id is known, then use the following command to show a plot with raw and windowed data,
replacing the string readID with the appropriate read id.
If you want to display a model curve as well, then you have to extract contig, start and end corresponding
to the read id and run the block above.

```shell
{ 
    # get raw data
    samtools view -b -e 'qname=="readID"' $mod_bam |\
      bedtools bamtobed -i stdin |\
      awk '{print $0 "\t" $4 "_" $1 "_" $2 "_" $3}' |\
      sed '1i\contig\tstart\tend\tread_id\tignore1\tignore2\talt_read_id' |\
      python get_raw_data_from_modBAM.py --piped-regions --alt-read-id-column $mod_bam |\
      awk '{print $1 " " $2 " " $2+1 " " $3 " rawDetect"}';
	  
    # window it
    samtools view -b -e 'qname=="readID"' $mod_bam |\
      bedtools bamtobed -i stdin |\
      awk '{print $0 "\t" $4 "_" $1 "_" $2 "_" $3}' |\
      sed '1i\contig\tstart\tend\tread_id\tignore1\tignore2\talt_read_id' |\
      python get_raw_data_from_modBAM.py --piped-regions --alt-read-id-column $mod_bam |\
      sed '1i\detectIndex\tposOnRef\tval' |\
      python get_mean_brdU_window.py --window 300 --thres 0.5 |\
      awk '{print $1 " " $3 " " $4 " " $2 " winDetect"}';
	          
} > plotting_and_short_analyses/plot_data

# plot data
cd plotting_and_short_analyses;
Rscript ./plot_one_read_w_win_or_model_if_needed.R plot_data plot.png
```

Note:
* One can choose to plot just raw data, just windowed data, or just model, or a combination by omitting
the appropriate lines above.
* Model parameters must already have been obtained to plot a model curve. Snippets above do not fit a model.
* To force a window boundary at a specific reference-genome position `num`, use the flag
`--forceWinBoundaryAtPos num` in the `get_mean_brdU_window.py` script.

Additional notes:
* If multiple model curves need to be plotted on the same plot, insert an NA after each model to ensure
the model curves get broken up. What I mean is: run the code to make the `plot_data` file below.
Then, run `tail -n 1 plot_data | awk '{print $1, $2, $3, "NA", $5}' >> plot_data`.
This will repeat the last model line, but replace the model value with NA.
Then add another model curve using the last python block below.
Then, add another NA line.
And so on and so forth.
* We are repeating a raw-data-obtaining command once.
There are ways to not repeat the command, but we are not going to bother.
If you are interested, look up how to use the `tee` command.

## Intermediate plotting workflows

If you are a somewhat experienced user and want to get finer control over plot elements, then use
the plotting scripts in this section. They are built on and coordinate the commands in the previous section.
Most or all of these plotting scripts are in the folder `plotting_and_short_analyses/`.

- To make rain plots from one read, use `plot_read.sh`.
- `plot_one_read_with_features.sh`, `plot_n_reads_given_detectIndices_and_sigmoids.sh` extract and process
data from model-fit files before feeding them to `plot_read.sh`.
- `convert_plots_to_latex.sh` gathers a set of rain plots into a pdf and displays information after each plot.
- `plot_forkSense_data_to_accompany_nascent_reads.sh` plots left and right fork probabilities per thymidine
on nascent reads.
 
### Plot read: make rain plots of one read with different features

In the most basic usage, the script takes a mod bam file, a read id, a contig, a start and an end as inputs,
and outputs a plot of the analogue modification probabilities per thymidine along the read in the specified region
and the associated data file to the specified output directory.

```bash
bash plot_read.sh sample.mod.bam readID chrII 1000 20000 output_dir
```

We add another layer now: a plot of windowed data.
Windows are of size 300 thymidines, and thymidines with a modification probability of 0.5
or more are considered modified, whereas ones with lower probabilities are considered unmodified.
A window boundary is forced to lie at the reference coordinate of 1100.

```bash
bash plot_read.sh sample.mod.bam readID chrII 1000 20000 output_dir 300 0.5 1100
```

We add yet another layer now: a plot of model predictions, and the fasta reference genome as an input.
Our plotted sigmoidal function has the following features/parameters: 
- it runs between the levels of 0.1 and 0.7 (these are asymptotic levels and may not be realized on the plot)
- the inflection point is at x = 10000 (inflection point is where the y value is (low + high)/2)
- the width of the sigmoid is 300 ( +- width from the inflection point is where the y value runs from
0.73 * low + 0.27 * high to 0.27 * low + 0.73 * high)
- the curve is plotted between x = 10000 and x = 20000
- the sigmoid runs from high to low with increasing reference coordinates. Flipped sigmoids use a negative width.

```bash
bash plot_read.sh sample.mod.bam readID chrII 1000 20000 output_dir 300 0.5 1100 fasta 0.7_0.1_10000_300_10000_20000
```

Read the comments in the script to understand how to plot multiple model curves per read, plot pauses,
plot annotations like fork directions, origins, terminations etc. using other command line arguments.  

### Plot read: make rain plots of one read, receiving pause information from a pause file

In the previous subsection, we saw how to plot many model curves per read using
many parameter strings in the format `high_low_offset_width_start_end`.
However, our model results are not stored in this format.
The script below converts model results stored in a pause file `$pause_file` to the parameter string
format and feeds it to `plot_read.sh`.
Consult comments in the script on what the pause file should look like
and what inputs such as `$feature_start` and `$feature_end` are.

```bash
bash plot_one_read_with_features.sh $read_id $mod_bam $feature_start $feature_end\
  $forksense_dir $pause_file $fasta_file $op_dir
```

### Make rain plots of read(s), receiving model fit information from a model fit file

When we fit variable width sigmoids to forks, we want to plot the data and the fitted sigmoids per read.
The script below converts the information in a model-fit file `$detect_index_sigmoid_file`
to the parameter string format and feeds it to `plot_read.sh`.
Consult the script on what the model-fit file should look like.

```bash
bash plot_n_reads_given_detectIndices_and_sigmoids.sh $detect_index_sigmoid_file $mod_bam $forksense_dir $fasta_file
  $mod_bam_left $mod_bam_right $op_dir
```

### Plot fork sense data to accompany plots of probability of analogue modification per read

An input folder contains many plots of data per read id.
The script locates these plots, extracts one read id per plot from the file name,
locates the forkSense data (raw left and right fork probabilities) for that read id if possible,
plots the forkSense data and sends it to the same folder.
We say 'if possible' in the sentence above as forkSense data may not be available for all reads.

```bash
sbatch plot_forkSense_data_to_accompany_nascent_reads.sh $modBAMForkSenseLFile $modBAMForkSenseRFile $opDir
```

* `$modBAMForkSenseLFile` mod bam file with left fork probabilities
* `$modBAMForkSenseRFile` mod bam file with right fork probabilities
* `$opDir` is the folder with several plots whose names are in the format `plot_<readID>.png`.

Can use `bash` in place of `sbatch` if need be.

### Collate several plots into pdf with associated information

The script gathers the following into one pdf:
- an input folder contains several plots identifiable by the read id
- an input data file contains rows of data, with each row corresponding to one read id
- a fork sense directory contains fork features per read id

```shell
bash convert_plots_to_latex.sh $modBAM_statistics $png_dir $forkSense_dir $other_plots_dir $prefix $feature_prefix
```

* `$modBAM_statistics` is a file with headers and tab-separated columns.
One header must be named `detectIndex` and this column must contain strings that contain the read id.
* `$png_dir` is any directory that contains many plots with names in the format `plot_<readID>.png` where readID is
any read id. The output pdf is sent to this directory.
* `$forkSense_dir` is the directory that contains left forks, right forks, origins etc. called by forkSense
* `$other_plots_dir` (optional parameter) is a directory that contains many plots with names
in the format `<prefix>_<readID>.png` where readID is any read id and prefix is as below.
* `$prefix` (optional parameter) is a prefix of plots in `$other_plots_dir` (see above).
* `$feature_prefix` (optional parameter) if fork features have already been collected in files
whose names are `$png_dir/<feature_prefix>_<readID>`, then those files are directly used,
instead of searching through files in `$forkSense_dir`.