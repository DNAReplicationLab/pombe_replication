# Example workflows using the codebase

## Running the main pipeline

Please consult `main_workflow.md` in this folder.
Also please read the section 'Preamble' below.

## Preamble

* Please ensure software tools (samtools, python etc.) have been installed, and that they contain
any necessary packages before running the snippets below.
* Some commands below may require a lot of memory or disk space depending on the input.
It is up to the user to provide these.
* Some of the commands here are bare-bones commands to illustrate basic usage of various scripts. Please write your own
scripts to coordinate these commands or use other already-written scripts in the repository that perform
such coordination.

## Convert detect to modBAM

BAM formats can include optional tags to store base-modification data. We call these modBAM files
in our notes and scripts. Convert detect to modBAM files using the snippet below. Only information
in the detect file (and associated fasta reference genome) is used to make the modBAM file.
Any other information such as cigar strings etc. from the bam files created earlier in the pipeline
will not be present in the modBAM file created below.

```shell
< $detectFile python DNAscentTools/convert_detect_to_modBAM.py --op "$modBAMPrefix".bam --tag T
```

## Convert modBAM to detect

```shell
samtools view -h $modBAM | python convert_modBAM_to_detect.py > outputDetectFile.detect
```

## Commands to extract raw data from modBAM

Each read id in a mod bam file is associated with modification probabilities at several bases.
We will discuss ways of extracting this data here.
If the mod bam file was made using our pipeline, we can use the `get_raw_data_from_modBAM.py` script
or the `modkit` program.
If the mod bam file was not made by us, then `modkit` works, but `get_raw_data_from_modBAM.py`
and other mod bam scripts in our repo will not, unless the custom `XR` and `XA` tags are added to the mod bam file
using the `add_XR_tag_to_modBAM.py` and `add_XA_tag_to_modBAM.py` scripts.
As all our scripts are based on the linux pipe `|` paradigm which allows for swaps of commands, 
we can easily switch between the two methods.

### Extract data from one read where read id, contig, start, end are known

Replace with suitable mod bam, read id, contig, start, end.

```shell
python get_raw_data_from_modBAM.py sample.mod.bam 00058fe1-e555-4a64-a41b-7f58fb7d6d6c chrXIV 243000 287000
```

NOTE: If a read id appears multiple times in the modBAM, then all or some or none
of these entries may be output.
If you run into this problem, please look at the logic in the `get_raw_data_from_modBAM` file.


### Extract data from many reads simultaneously where read id, contig, start, end per read are known

Create a file with space- or tab-separated columns with headers read id, contig, start, end in any order.

```shell
< regions_file python get_raw_data_from_modBAM.py --piped-regions sample.mod.bam
```

### Extract data from one or many reads where only read ids are known

We use samtools and bedtools to fill in contig, start, end per read. Prepare a file whose name is
`$readIDs` in the snippet below which contains a one column of read ids with no header.

```shell
samtools view -b -N $readIDs $modBAM |\
    bedtools bamtobed -i stdin |\
    awk '{print $0 "\t" $4 "_" $1 "_" $2 "_" $3}' |\
    sed '1i\contig\tstart\tend\tread_id\tignore1\tignore2\talt_read_id' |\
    python get_raw_data_from_modBAM.py --piped-regions --alt-read-id-column $modBAM
```

### Extract data from all reads in the mod bam file

From the previous subsection, just remove the `-N $readIDs` from the samtools command.
Alternatively, combine the samtools and bedtools command into `bedtools bamtobed -i $modBAM`,
which can be fed into the awk command.

### Extract data from one read where only read id is known

The method used for many reads can be used here as well.
Alternatively, use the following, replacing readID with the actual readID.
NOTE: it is assumed that there is only one entry per readID.

```shell
samtools view -b -e 'qname=="readID"' $modBAM |\
    bedtools bamtobed -i stdin |\
    awk '{print $0 "\t" $4 "_" $1 "_" $2 "_" $3}' |\
    sed '1i\contig\tstart\tend\tread_id\tignore1\tignore2\talt_read_id' |\
    python get_raw_data_from_modBAM.py --piped-regions --alt-read-id-column $modBAM
```

### Extract data from many forks 

Our notation for forks is `readID_contig_start_end_orn_dirn_startFork_endFork`.
`readID, contig, start, end, orn` are pieces of information about the alignment on which the fork was called.
It must be clear by now what the first four quantities mean.
`orn` must be `fwd` or `rev`.
`dirn` which means fork direction and must be `L` or `R`.
`startFork` and `endFork` are locations on the reference genome corresponding to a fork call.
Our convention is that `startFork < endFork` irrespective of `dirn`.

Create a file with one column and header `alt_read_id`, consisting of one fork per row in the notation above.
The command below extracts raw data between the fork coordinates.

```shell
< fork_file python get_raw_data_from_modBAM.py --piped-regions --fork-info-in-alt-read-id sample.mod.bam
```

## Other commands to process detect data in modBAM form

### Output mean analogue density per read averaged over entire read

We are appropriating a script used for the different purpose of identifying nascent reads,
by setting the nascent threshold to zero and averaging over the entire read (using a window size of 0).

```shell
samtools view $modBAM |\
    python get_modBAM_nascent_reads.py --window 0 --nascentThres 0 --showVal 
```

Another method is to use scripts mentioned before to extract data from all reads of the mod bam file,
and then use the windowing script with the default window size of NaN which will average over the entire read.

```shell
bedtools bamtobed -i $modBAM |\
  awk '{print $0 "\t" $4 "_" $1 "_" $2 "_" $3}' |\
  sed '1i\contig\tstart\tend\tread_id\tignore1\tignore2\talt_read_id' |\
  python get_raw_data_from_modBAM.py --piped-regions --alt-read-id-column $modBAM |\
  sed '1idetectIndex\tposOnRef\tval' |\
  python get_mean_brdU_window.py --thres 0.5
```

### Output mean analogue density per many non-overlapping windows of a fixed number of thymidines per read

Use the command above, but set a window size in the last command using the flag `--window`.
For e.g. `python get_mean_brdU_window.py --thres 0.5 --window 300`.

### Output mean analogue density per sliding windows of a fixed number of thymidines per read

Now, we want to set up sliding windows of a fixed number of thymidines, and a step size of 1 thymidine. 
Use the command above, but set an additional flag `--sliding`.
For e.g. `python get_mean_brdU_window.py --thres 0.5 --window 300 --sliding`.

### Window detect data corresponding to one read at one genomic location

```shell
python get_raw_data_from_modBAM.py $modBAM $readID chrXIV 243820 287092 |\
    sed '1idetectIndex\tposOnRef\tval' |\
    python get_mean_brdU_window.py --window 300 --thres 0.5
```

* To set window boundary at right end and report data from last window to first window,
use the flag '--rev'.
* Use the flag '--infer' to automatically infer if data is to be
reported first to last window or vice versa depending on whether the strings "\_R\_" or "\_L\_"
are present in detectIndex.

### Extract and average detect data corresponding to all reads at one genomic location

The following example uses the window on chrX between 10000 and 11000 as an example.
It extracts all reads which pass through all the interval chrX:10000-11000 and averages
modification density over this interval per read and reports one value per read.
NOTE: Any window size can be used. It is a coincidence that the window size here is 1 kb.

```bash
contig=chrX
start=10000
end=11000

samtools view -e "pos<=$start&&pos+rlen>=$end" $modBAM $contig:$start-$end |\
    awk '{print "'"$contig"'",$1,'"$start"','"$end"'}' |\
    sed '1i\contig read_id start end' |\
    python get_raw_data_from_modBAM.py --piped-regions $modBAM |\
    sed '1idetectIndex\tposOnRef\tval' |\
    python get_mean_brdU_window.py --thres 0.5
```

Sample few lines from output (not real data)

```text
f034b4d7-891e-4c77-b19b-ed86593f0833 0.064 10002 10998
109e3a1c-e95d-4d1a-a0fd-37a276fe333e 0.047 10002 10998
```

### Get mean of detect data over region without windowing and with or without thresholding

```shell
python get_raw_data_from_modBAM.py $modBAM $readID chrI 0 2000 |\
    sed '1idetectIndex\tposOnRef\tval' |\
    python get_mean_brdU_window.py --thres 0.5
    
python get_raw_data_from_modBAM.py $modBAM $readID chrI 0 2000 |\
    sed '1idetectIndex\tposOnRef\tval' |\
    python get_mean_brdU_window.py  
```

### Extract detect headers from modBAM files

Extract all headers

```shell
samtools view $modBAM | awk -F 'XA:Z:' '{print $2}' |\
    awk '{print $1}' | sed 's/_/ /g' | sed 's/^/>/g'
```

Extract header corresponding to one read id

```shell
samtools view -e 'qname=="read-id-here"'  $modBAM |\
    awk -F 'XA:Z:' '{print $2}' |\
    awk '{print $1}' | sed 's/_/ /g' | sed 's/^/>/g'
```

### Extract nascent read ids from modBAM files

Strategy 1: if a read has at least one window of size 300 thymidines with analogue density >= 0.05 after
thresholding, then that read is regarded as nascent. Windows are not sliding, they are non-overlapping.

```shell
samtools view $modBAM | python get_modBAM_nascent_reads.py --window 300 --thres 0.5 --nascentThres 0.05 --tag T
```

Strategy 2: if overall analogue density in the read is >= 0.05 after thresholding, then that read is regarded
as nascent. Note that in this script, a window size of 0 means average over the entire read.

```shell
samtools view $modBAM | python get_modBAM_nascent_reads.py --window 0 --thres 0.5 --nascentThres 0.05 --tag T
```

## Working with modBAM files made by other workflows

Other research projects may not use DNAscent to generate a .mod.bam file.
We need to make some changes to such a file before our tools can be used on it.
If the file is not reference-anchored, then expect a lot of our tools to not work on it.
This may be fixed in the future.
If the file is reference-anchored, then execute the following commands to add XR and XA tags to them,
which a lot of our tools require.

```shell
samtools view sample.mod.bam | python add_XR_tag_to_modBAM.py |\
  python add_XA_tag_to_modBAM.py | samtools view -b > sample.output.mod.bam
```

The XR tag's value per read is an integer obtained by converting the first seven characters of a read id
from hexadecimal to decimal.
The XA tag's value per read is a string `readID_contig_start_end_orn`. contig, start, and end refer
to positions on the reference genome whereas orn is fwd/rev.

## Fit pauses to data from one read

Run the cut-and-align procedure to find pauses on thresholded (at 0.5) raw data from a read
using a sigmoid that runs from 0.7 to 0.1 with increasing x and width 3 kb as a reference right-moving
fork.

NOTE: There are two options to specify fork direction. (1) A positive (negative) width parameter
means a right-moving (left-moving) fork. (2) The detect index must contain the string `_L_` or `_R_`.

```shell
python get_raw_data_from_modBAM.py $bamfile $readid $contig $start $end |\
    sed '1i\detectIndex\tposOnRef\tprobBrdU' |\
    python get_model_pauses_raw_sgm.py --thres 0.5 --width 3000 --low 0.1 --high 0.7 --iter 20 --tol 0.03
```

Command above produces five output columns: detectIndex, pauseDuration, pauseSite, gof, leftCutPt.
To obtain best-fit pause site and assess if the site and duration are such that we are confident
of pause detection, redirect the output from above to the two additional scripts below as shown,
noting how the first get raw data command has been slightly modified.

```shell
  python get_raw_data_from_modBAM.py $bamfile $readid $contig $start $end False ${readid}_${contig}_R_|\
    sed '1i\detectIndex\tposOnRef\tprobBrdU' |\
    python get_model_pauses_raw_sgm.py --thres 0.5 --width 3000 --low 0.1 --high 0.7 --iter 20 --tol 0.03
    python process_model_pause_results_raw_sgm.py --width 3000 --low 0.1 --high 0.7 |\
    python process_model_pause_results_analysis.py \
        --modbam "$bamfile"\
        --modbam_left "$modbam_left" --modbam_right "$modbam_right"\
        --keep_cols_requested keep_exp_gof_v2,keep_fl,keep_fp,keep_ep,keep_lp,keep_fs
```

Commands above produce many output columns. Optional: to improve readability, you can use the package miller (`mlr`)
to convert CSV into JSON.

## Fit a pause per read over many reads simultaneously

Prepare a file with space or tab-separated columns with the headers: read_id contig start end.
Headers can be in any order and additional headers are ok.

Execute the pause-finding commands as above, with the raw data retrieval step as below.

```shell
< $regions_file python get_raw_data_from_modBAM.py --piped-regions $bamfile
```

## Fit one sigmoid to collection of windowed data from many reads simultaneously

Now, we want to fit sigmoids to windowed analogue probabilities than their raw counterparts.
Prepare a file regions_file with the columns contig start end read_id alt_read_id where alt_read_id
is a unique label per fork which can be anything as long as it contains the strings `_L_` or
`_R_` to denote the fork is leftward or rightward. We've assumed the file has no headers and is
space-separated in the snippet below.

* Set a variable paramsFitFile where best-fit parameters etc. can be reported.
* Set variable modbam to location of modbam file
* One of the window sizes below is in thymidines whereas the other is in bases.
* There are outputs to stdout as well.

```shell
< $regions_file sed '1i\contig start end read_id alt_read_id' |\
    python get_raw_data_from_modBAM.py --piped-regions --alt-read-id-column $modbam |\
    sed '1idetectIndex\tposOnRef\tval' |\
    python get_mean_brdU_window.py --window 300 --thres 0.5 --infer |\
    sed '1i\detectIndex mean_brdU start end' |\
    python get_model_fit.py --win 1000 --fitFile $paramsFitFile --xMax 45 --winTol 0.1
```

## Fit one sigmoid to collection of windowed data from many forks

Use the script `run_get_model_fit.sh`, which runs logic similar to the lines in the previous section.
Inputs are left and right forks listed in the format used by forkSense, and a mod bam file.
Outputs are sent to the specified location.
Parameters can be adjusted in the script.

## Prepare a regions file from forkSense fork direction calls

### Preamble

A sample fork sense fork direction file looks like

```
chrI 8976 13925 0bd9a5a0-1f8f-49ca-be79-de597f6c32fb chrI 0 22438 fwd
chrI 165 4939 750f5072-9728-4daf-a09b-e98cabf5dadf chrI 0 5104 fwd
```

The information about whether these are left- or right- moving forks is not
stored in the file but may have been recorded in the filename or elsewhere.
NOTE: fork direction cannot be inferred from the second and third columns, as
the second column is always less than the third column.

### Convert into a `regions_file` for usage in `get_raw_data_from_modBAM.py`

Shell variables below corresponding to different filenames must be set suitably

```shell
{
< $fork_sense_left_file awk 'BEGIN {OFS=" "}{print $1, $2, $3, $4, sprintf("%s_%s_%s_%s_%s_%s_%s_%s",$4,$1,$6,$7,$8,"L",$2,$3)}' 
< $fork_sense_right_file awk 'BEGIN {OFS=" "}{print $1, $2, $3, $4, sprintf("%s_%s_%s_%s_%s_%s_%s_%s",$4,$1,$6,$7,$8,"R",$2,$3)}' 
} > $regions_file
```

## Measure an externally-supplied quantity vs analogue density in genomic windows

Goal: we have a dataset obtained from elsewhere containing one measurement per genomic window.
We want to measure mean analogue density in each genomic window across all reads so that
we can do further analyses elsewhere on the relationship between the prior measurement and the
analogue density. For the purposes of this section, let us assume the measurement made
elsewhere is the median replication time. The data here is for illustrative purposes only;
do not assume that it is real data.

Scripts used

```text
plotting_and_short_analyses/plot_trep_vs_BrdU_mean.R 
get_agg_brdU_ref_coord_int.py
run_get_agg_brdU_ref_coord_int.py
```

The logic in the scripts is explained henceforth.
Our data file obtained from elsewhere looks like

```bash
head $dataFile
```

```text
chrI    3000    3999    +       50.9200
chrI    3999    4999    +       51.7200
chrI    4999    6000    +       53.2800
chrI    6000    6999    +       53.3600
chrI    6999    7999    +       53.4800
chrI    7999    8999    +       53.0700
chrI    8999    9999    +       52.5100
chrI    9999    10999   +       52.2000
chrI    10999   11999   +       52.9000
chrI    12999   13999   +       50.7300
```

Here, we obtain mean analogue density and total number of modified and unmodified bases per each genomic window
from a .mod.bam file that was made using the tools in this repository.

```bash
< $dataFile sed '1i\contig\tstart\tend\torn\ttrep' |\
    python get_agg_brdU_ref_coord_int.py --modBAM $modBAM --thres 0.5
```
 
```text
chrI    3000    3999    +       50.9200 0.0620  8734
chrI    3999    4999    +       51.7200 0.0592  4710
chrI    4999    6000    +       53.2800 0.0545  1732
chrI    6000    6999    +       53.3600 0.0560  9087
chrI    6999    7999    +       53.4800 0.0620  6381
chrI    7999    8999    +       53.0700 0.0609  3124
chrI    8999    9999    +       52.5100 0.0614  5678
chrI    9999    10999   +       52.2000 0.0624  1021
chrI    10999   11999   +       52.9000 0.0596  8129
chrI    12999   13999   +       50.7300 0.0675  2683
```

If this output is sent to a file, it can be visualized using the R script mentioned earlier in the section.

## Measure average analogue density per regular genomic window across forks

Goal: Tile the genome into regular windows and measure average analogue density within the window across
all forks that pass through 100% of the window.

### Program logic

We tile the genome into regular windows of a given size.
We average every fork along its direction of movement in *thymidine* window sizes 3 times the size of the
regular window, which corresponds to genomic window sizes of approximately 10 times the size of the regular window.
So, we have two sets of data now:
- a set of windows w_i from the genome where i is an index along the genome.
- a set of windows and analogue average densities a_ij, b_ij from the forks where i is an index that runs along the
fork, and j is an index that labels forks.

For each w_i, we calculate mean(b_kl) for all k and l where w_i is completely contained in the window a_kl.

We need to think carefully about the window size parameter that is input to the program.
If you set the parameter to say a 100 bp, that means the genome is tiled in 100 bp windows and the forks are averaged
using approximately 1000 bp windows. You have to think about your specific question and set the window size
accordingly.


### Sample output

A few lines from the output may look like:

```text
chrI 1600 1700 0.283 13
chrI 1700 1800 0.2592 15
```

We can visualize the output using IGV as it is in the bedgraph format (only the first four columns are used by IGV).
The first line means that 13 forks passed through the window chrI:1600-1700 and their mean analogue density is 0.283.

### Command

```bash
sbatch -o $output_file calculate_brdu_amount_vs_genome_window_on_forks.sh $fasta_index_file $window_size $mod_bam_file \
         $left_fork_file $right_fork_file
```

Let's look at what the script does through an example.
NOTE: \ means the command continues on the next line.

### Example

```bash
sbatch -o /path/to/output/file.bedgraph calculate_brdu_amount_vs_genome_window_on_forks.sh /path/to/sacCer3.fa.fai \
  100 /path/to/sample.mod.bam /path/to/leftForks_DNAscent_forkSense.bed /path/to/rightForks_DNAscent_forkSense.bed
```

We use
- analogue modification data from `/path/to/sample.mod.bam`. The file must have been sorted and indexed.
- forkSense-called left and right fork data from `/path/to/leftForks_DNAscent_forkSense.bed` and
`/path/to/rightForks_DNAscent_forkSense.bed`.
- A fasta index file of the reference genome `/path/to/sacCer3.fa.fai`. If you have a fasta file, you can generate
the index file using `samtools faidx`.

We tile the genome into windows of size set by the corresponding input parameter, which is 100 bp in this example.
Read the program logic section to understand how to choose a good window size.
Output is written to the bedgraph `/path/to/output/file.bedgraph`.

## Measure average analogue density per regular genomic window across features

Please read the previous section entitled 'Measure average analogue density per regular genomic window across forks'
first. Here, we perform the same calculation but for features other than forks called by forksense (e.g. origins,
terminations). The command is

```bash
sbatch -o $output_file calculate_brdu_amount_vs_genome_window_on_features.sh $fasta_index_file $window_size\
   $mod_bam_file $feature_file
```

## Measure fork/feature count per regular genomic window

Goal: Tile the genome into regular windows and count all forks per window that pass through 100% of the window.
We can also count features other than forks called by forksense by using suitable input files (e.g. origins,
terminations).

### Program logic

We tile the genome into regular windows of a given size. So, we end up with:
- a set of windows w_i from the genome where i is an index along the genome.
- a set of forks with start s_i and end e_i, where i is an index that that labels forks.

For each w_i, we count forks k where w_i is completely contained in the window s_k to e_k.
You have to think about your specific question and set the window size for fork counts accordingly.

### Sample output

A few lines from the output may look like:

```text
chrI 1600 1700 13
chrI 1700 1800 15
```

We can visualize the output using IGV as it is in the bedgraph format.
The first line means that 13 forks passed through the window chrI:1600-1700.

### Command

The script may not be in the root directory of the repository, but might be in some other directory like
`plotting_and_short_analyses`. 

```bash
bash calculate_feature_count_vs_genome_window.sh $fasta_index_file $window_size $feature_file_1 $feature_file_2 ... \
     > output_file.bedgraph
```

Let's look at what the script does through an example.

NOTE: 
- \ means the command continues on the next line,
- ... means that there can be more than one feature file,
- \> means that the output which is normally sent to the terminal will be sent to the file `output_file.bedgraph`,
creating the file if it does not exist and overwriting it if it does exist.

### Example

```bash
bash calculate_feature_count_vs_genome_window.sh /path/to/sacCer3.fa.fai 100 \
  /path/to/leftForks_DNAscent_forkSense.bed /path/to/rightForks_DNAscent_forkSense.bed > /path/to/output/file.bedgraph
```

We use
- forkSense-called left and right fork data from `/path/to/leftForks_DNAscent_forkSense.bed` and
`/path/to/rightForks_DNAscent_forkSense.bed`.
- A fasta index file of the reference genome `/path/to/sacCer3.fa.fai`. If you have a fasta file, you can generate
the index file using `samtools faidx`.

We tile the genome into windows of size set by the corresponding input parameter, which is 100 bp in this example.
Read the program logic section to understand how to choose a good window size.
Output is written to the bedgraph `/path/to/output/file.bedgraph`.


## Query fork features associated with a read id

Our pipeline calls fork features such as origins, terminations, fork directions etc. per read id.
If the pause detection part of the codebase is run, additionally pauses can be detected per read id.
If one wishes to query the features associated with a read id, the procedure is as follows.

Case I: Pause detection not performed, but fork sense etc. have been run.
Use the following command.

```bash
bash get_feature_information.sh $read_id /dev/null $forkSenseDirectory
```

`$forkSenseDirectory` is where several '.bed' files were deposited by forkSense.
`/dev/null` is a special, empty file in a linux system.
We use such a file in place of a file containing detected pauses.

Sample output (for illustrative purposes only. Do not assume this is real data)

```text
1502846 1503345 origin
1473702 1473927 origin
1519112 1522710 leftFork
1490807 1502846 leftFork
1468214 1473702 leftFork
1473927 1488934 rightFork
1503345 1517847 rightFork
1517448 1519402 termination
1488467 1491448 termination
```

Case II: Pause detection performed in addition to running fork sense etc.
Use the following command.

```bash
bash get_feature_information.sh $read_id $pause_file $forkSenseDirectory
```

`$pause_file` is tab-separated with headers.
The columns `detectIndex`, `pauseSite`, and `keep*` are used in the script.
`detectIndex` is any string that contains the read id in it,
`pauseSite` is self-explanatory,
`keep*` are columns that are set to 'True' if the pause site is expected to be real.
There are many such columns (as implied by the wildcard '*') as several
independent methods may be used to determine if a pause site is real.

Line(s) like the following may appear at the end of the output.

```text
1470463 1470463 pause
```

## Filter forksense origin and fork files based on length and contig criteria

Our goal here is to create subsets of forksense origin and fork files based on criteria of fork length,
alignment length, and contig.

### Input files

Forksense creates files called x_DNAscent_forkSense.bed where x = origins, leftForks, rightForks, or terminations.
Note that names may be different, and that output files depend on the options and the version of forkSense used.
If names are different, then the script may not work.

These files look like the following (don't assume this is real data).

```text
chrI 1155 9507 a437199d-5381-49ac-bd0e-d8283255f1aa chrI 0 10181 fwd
chrI 291 4930 924c09bb-fc9a-4169-9156-ffe8630bfe4a chrI 0 5291 fwd
```

The above lines mean a feature exists between the positions in column 2 and column 3 on contig 'chrI' on the molecules
of the respective read ids.

### Command

We want to preserve only those origins and forks that lie on some contigs, and that are on forks and/or alignments
of a certain length. The command to be used is as follows. Note that the script may not be in the root directory
of the repository and may be in a subdirectory like `specialized_analyses`.

```bash
bash filter_forksense_files_on_length_and_noChrM.sh $forksense_dir $min_fork_length $min_align_length $output_dir $chrM_name
```

Let's look at what the command does through some examples

### Example 1

```bash
bash filter_forksense_files_on_length_and_noChrM.sh /path/to/forksense/dir 5 10 /path/to/output/dir chrI
```

Here, the forksense files like origins.bed etc. are present in the directory `/path/to/forksense/dir`.
We want to preserve only those features that are associated with forks and alignments of length at least 5 kb and 10 kb
respectively.
We also want to exclude features on contig 'chrI'.
Output files are sent to the directory `/path/to/output/dir` and are called
x_DNAscent_forkSense.forkLen.5kb.alignLen.10kb.noChrI.bed where x = origins, leftForks, or rightForks.
Note that we do not filter terminations.

### Example 2

```bash
bash filter_forksense_files_on_length_and_noChrM.sh /path/to/forksense/dir
```

Like in the previous example, the forksense files are present in the directory `/path/to/forksense/dir`.
The other arguments are not provided, so the script uses default values.
Please see the script for these default values.

## Convert a bed file into coverage information

Goal: Tile the genome into windows of a given size and count the number of overlapping bed
entries per window. We only count a bed entry if the overlap is 100% with a genomic window.

### Program logic

We tile the genome into windows of size set by the corresponding input parameter.
For each genome interval, we count the number of bed intervals whose start and end positions are such that the
bed interval passes through 100% of the genome interval.

### Command and example

```bash
bash convert_bed_to_coverage.sh /path/to/genome.fa.fai 100 /path/to/file.bed > /path/to/output/file.bedgraph
```

We use the following input files:
- /path/to/genome.fa.fai: a tab-separated file with many columns and no header. We'll use only the first two columns of
contig name and contig length. These files are usually found in the same directory as the genome fasta file.
If not available, then use the command `samtools faidx /path/to/genome.fa` to generate the file.
- We use a window size depending on the application/context. Here, we use 100 bp as an input parameter.
- /path/to/file.bed: a tab-separated file with no column names. We'll use only the first three columns.
- Normally, output is sent to stdout. We redirect it to a file `/path/to/output/file.bedgraph` using the `>` operator.
This file is created if it does not exist, and overwritten if it does exist.

#### Invert option

Optionally, we can invert the overlap criterion i.e. instead of requiring the bed entry to overlap 100% with the genome
window, we can require the genome window to overlap 100% with the bed entry. Such a command looks like

```bash
bash convert_bed_to_coverage.sh /path/to/genome.fa.fai 100 /path/to/file.bed --invert > /path/to/output/file.bedgraph
```

We can use this command if the intervals in the bed file are much smaller than 100 bp and our context/application
requires such an analysis.

### Sample output

A few lines from the output may look like:

```text
chrI 1600 1700 13
chrI 1700 1800 15
```

We can visualize the output using IGV as it is in the bedgraph format.
The first line means that 13 bed intervals passed through the window chrI:1600-1700.

## DNAscent qc

See `dnascent_qc.md` for details on how to run the DNAscent qc script.