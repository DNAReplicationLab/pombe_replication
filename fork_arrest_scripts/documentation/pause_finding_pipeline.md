# Non-repetitive genome

You can run `bash get_sgm_dnascent_pauses.sh` to find single-molecule pauses.
You must have already obtained analogue modification calls, replication track
calls, and fit a reference sigmoid to run this step;
consult `main_workflow.md` if you do not have these files.
You need to prepare an input json file like the following excerpt.
Low, high, and width refer to the parameters of the reference sigmoid.

NOTE: we generally use cutoffs on fork length and alignment length on the forks
we feed to the pause-finding pipeline; a typical threshold is 3 kb and 30 kb respectively.
You can use the program `filter_forksense_files_on_length_and_noChrM.sh` to perform this filtration.

```json
{
   "left_forks": "/path/to/left/fork/calls.bed",
   "right_forks": "/path/to/right/fork/calls.bed",
   "mod_bam": "/path/to/mod.bam",
   "output_dir": "/path/to/output/dir",
   "output_file_prefix": "pauses_08may24",
   "low": 0.11,
   "high": 0.65,
   "width": 3060,
   "mod_bam_left": "/path/to/forkSense/left/mod.bam",
   "mod_bam_right": "/path/to/forkSense/right/mod.bam",
   "pipeline_go_get_model_pauses": 1,
   "pipeline_go_perform_further_analysis": 1
}
```

# rDNA

For finding pauses in the rDNA, we use a simpler pipeline. Please run the script `run_rDNA_detectSummary.sh`
with the appropriate input files.