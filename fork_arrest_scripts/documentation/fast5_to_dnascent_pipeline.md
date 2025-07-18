Here we describe how to go from nanopore currents to DNAscent modification calls and
replication track annotations using forkSense. If you already have modBAM files
(that you have generated earlier or downloaded from somewhere), you probably do
not need this step.

The main script you run is called `pipeline_fast5_pod5_to_dnascent_dorado.sh`.
You need an input json file that looks like the following excerpt.
The inputs below are all example inputs; you must fill them with locations to real files or realistic parameters.
You can see the comments in that file to see what inputs it needs in detail.

The scenario shown below is that you have a dataset where basecalling was performed on the device.
You want to use the fast5 and fastq files generated on device as input to the pipeline.
You also need a linear reference genome, and an HDF5 plugin needed by DNAscent v2.
If you have a slightly different scenario e.g. you want to basecall again, you
have to set the pipeline flags below differently.

```json
{
  "ref_genome_fasta": "/path/to/fasta/file.fa",
  "min_len": 1000,
  "min_qual": 20,
  "only_prim": 1,
  "dataset_prefix": "20210315_EXP_ONT_XY_abc123",
  "fast5_dir": "/path/to/dir/w/fast5/files",
  "fastq_dir": "/path/to/dir/w/fastq/files/fastq_*",
  "fastq_file": "/path/to/dir/w/fastq/files/fastq_*/*.fastq*",
  "guppy_sequencing_summary_file": "/path/to/sequencing_summary.txt",
  "hdf5_plugin_path": "/path/to/hdf5/plugin/needed/by/dnascent",
  "pipeline_go_align": 1,
  "pipeline_go_filter_sort_before_dnascent": 1,
  "pipeline_go_pycoQC": 1,
  "pipeline_go_dnascent_index": 1,
  "pipeline_go_detect": 1,
  "pipeline_go_make_modbam": 1,
  "pipeline_go_forkSense": 1,
  "pipeline_go_dnascent_qc": 1,
  "pipeline_go_make_forkSense_modbam": 1,
  "pipeline_go_plot_random_nascent_reads": 1,
  "pipeline_go_plot_random_nascent_reads_forkSense": 1,
  "pipeline_go_plot_random_nascent_reads_collate_pdf": 1,
  "pipeline_go_plot_filtered_random_nascent_reads": 1
}
```

The pipeline produces many outputs, such as:
- a mod bam file with modification calls per read
- many forkSense files that call replication tracks on these molecules
- a dnascent qc pdf, pycoQC analysis; these contain many useful plots like
read length distribution, analogue density distribution etc.; see `dnascent_qc.md`.

NOTE: If you are having a look at the data from the manuscript, it is worth mentioning that 
the basecalled sequences and modification-called modBAM files have been deposited in suitable
archives. So, you do not need to run this pipeline script again.
In case you want to run the basecalling and modification calling steps again,
you will need a proprietary basecaller, guppy, and the
older version v2 of the modification-calling software DNAscent.