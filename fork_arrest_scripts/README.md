# Fork arrest pipeline

Pipeline of the fork_arrest project, where we go from nanopore data
-> DNA sequence with BrdU substituted for dT in some positions, and perform analyses on these data
such as finding replication fork pauses and correlating their locations with other datasets.
This repository is primarily associated with a manuscript in preparation (Thiyagarajan et al.) and
is also used by another manuscript in preparation from the same research group (Díez Santos et al.).

## Preamble: dataset

You may want to run the scripts here on your data, or you may want to repeat analyses from the paper.
If it is the latter, then please download the dataset corresponding to the paper and have a look there first.

## Preamble: software and preparation

NOTE: If you are from the Earlham Institute and run programs on its HPC, then this section does not apply to you.

### Software

Please download the repository and any associated submodules such as DNAscentTools.
You will probably need to execute commands similar to the ones below.

```bash
# Exercise
git clone <url_to_this_repo>
# NOTE: if the above doesn't work, it is worth trying adding .git to the end of
# the URL i.e.
# git clone <url_to_this_repo>.git
cd <name_of_the_repo>
git submodule update --init --recursive # This should get associated repos like DNAscentTools 
```

The scripts in this repository uses a lot of other software packages.
We provide a singularity definition file here `singularity.def` that can be used to build a singularity container
with most of the software needed to run the scripts in this repository.
We had to omit some due to reasons mentioned in the header of the singularity definition file.
Or you can make a list of the packages from that file and install them on your own or check if you already have them.

### Preparation: load_package.sh

You've to edit `load_package.sh` suitably to point to the singularity image you've installed above.
For example, let's say you've installed the singularity image at `/a/b/c/singularity_image.img`
then the R lines in the script should be edited to read as follows:
```bash
    elif [ "$var" == "-R" ]
    then
        # R
        R() {
            singularity exec /a/b/c/singularity_image.img R --no-save "$@";
        }
        export -f R;
        Rscript() {
            singularity exec /a/b/c/singularity_image.img Rscript "$@";
        }
        export -f Rscript;
```
You've to do this for all the software packages that you use in this repository.

If you do not use singularity and have installed packages directly on your system, then you can just delete most of
the text in the file. This will render statements like `source load_package.sh -python` non-operational,
but that is fine as long as you have the required packages installed on your system.

### Preparation: config.sh

You must make a file `config.sh` with the following contents:

```bash
#!/bin/bash
declare -A config
config["dataDir"]="/where/you/have/downloaded/the/papers/dataset/or/made/your/own/data"
config["scratchDir"]="/a/temporary/directory"
config["name"]="Your name"
config["email"]="Your email"
config["job_output_logs"]="/path/to/job/output/logs"
```

## Introduction

There are three types of scripts in this repository:
- pipeline scripts that run the main analysis steps. These are described in the `Running pipelines` section below.
- scripts that are called by the pipeline scripts to perform specific tasks.
- other scripts that are used to perform specific analyses or to generate plots.

All the scripts in the repository have simple names and comments at the top that describe what they do.
e.g. `bed_region_pause_report.sh` generates a report of replication fork pauses in a region of the genome
specified by a bed file.

## Running pipelines

There are three main pipelines in this repository:
- The command `bash pipeline_fast5_pod5_to_dnascent_dorado.sh` runs the pipeline that converts
  raw nanopore data into modification calls.
- The script `subsample_forks_fit_cg_model_and_plot.sh` fits sigmoids to a subsampled set of
  forks. You can run this a few times to get a best-fit set of parameters per run, and then
  average these to obtain the reference sigmoid parameters you will use in the next step.
- `bash get_sgm_dnascent_pauses.sh` picks up from the previous step and finds single molecule
  replication fork pauses.
- `bash pipeline_sgm_pauses_winnow_and_analyse.sh` picks up from the previous step and performs
  further analyses on the replication fork pauses and correlates them with other datasets.

The three pipeline scripts are designed to launch jobs in a slurm-based HPC environment.
Scripts can be retooled to run in a local computer if necessary.
Each pipeline runs several steps each of which can be turned on or off.

The pipeline scripts require several input parameters, detailing where different input files are,
and what options to use for each step. These input parameters are largely specified in the json format.
The format of such files is described in the comments of the corresponding scripts.

Some more information about running the pipelines is in the `documentation/` directory of this
repository e.g. `pause_analysis_pipeline.md`.

## File formats

The main file formats used in the scripts are fast5, pod5, fastq, bam, modbam, bed, fasta, detect, forkSense, pauseFile
in addition to standard formats such as pdf, png etc. These formats are all well-known in the genomics field and
please google them if you are not familiar with them. The unfamiliar formats are: detect and forkSense which are
used by DNAscent detect and pauseFile which is a tab-separated value format used by us. Please refer to the
[documentation](https://dnascent.readthedocs.io/en/latest/detect.html) of the DNAscent repository for more
information on detect and forkSense. The pauseFile format is described in the `documentation/` directory.

## Using `sbatch` to run cluster commands

To run shell scripts on an HPC, use the sbatch command e.g.:
`sbatch limited_plot_read.sh $mod_bam $read_id $window_size $output_dir`.
One can add more flags to customize the job.
For example: `sbatch --mail-type=END,FAIL --mail-user=you@me.com limited_plot_read.sh ...`
will email "you@me.com" whenever the job ends in success or failure.
(the string "..." means other parameters are present).
Consult a guide on how to use slurm for more details.