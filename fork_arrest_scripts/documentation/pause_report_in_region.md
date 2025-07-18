Given a list of pauses and a list of regions on the genome, we want to answer
questions like 'how many pauses occur over this region?', 'how many are
head-on and how many are co-directional pauses?', 'how many pauses do
we expect in this region if pauses were uniformly distributed?',
'can I visualize the single-molecule data behind these pauses?' etc.
The script `bed_region_pause_report.sh` does just that.

This script produces a pdf file with a bunch of
statistics on whether pauses are enriched at a particular region in a bed file.
The script can optionally produce per-read plots of reads that overlap with
the region. You can find in this folder a sample report `sample_pause_report.pdf`
made from a bed file containing just one region (please do not assume this file
was made from real data).

The script takes a json file as input which looks like the
following (all are dummy values below, please fill in appropriate values):

```json
{
   "mod_bam": "/some/file.detect.mod.sorted.bam",
   "forksense_dir": "/some/directory",
   "pause_file": "/some/pause/file",
   "bed_file": "/some/file.bed",
   "op_dir": "/some/directory",
   "dataset": "some_dataset",
   "feature": "some_feature",
   "division": "0",
   "analysis_label": "some_label",
   "mod_bam_left": "/some/file.forkSense.mod.left.sorted.bam",
   "mod_bam_right": "/some/file.forkSense.mod.right.sorted.bam",
   "alignment_file": "/some/file.sorted.bam",
   "fasta": "/some/file.fna",
   "n": 20,
   "prefix": "collection",
   "prefix_plot_option": false,
   "relative_direction_option": "all",
   "delete_new_mod_bam_made_on_the_fly": false,
   "genome_size_bp_optional": 1000000,
   "AT_ratio_optional": 0.50,
   "pause_sensitivity_bedgraphs_prefix_optional": "/some/file.pause_sensitivity"
}
```

The regions of interest are in `/some/file.bed`.
The file has to be at least a six-column bed file (this may change in the future).
For more information about the parameters, please see the comments in the script.

If you do not have the necessary input files, please consult `main_workflow.md`.