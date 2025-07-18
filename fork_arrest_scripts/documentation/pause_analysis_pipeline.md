Here we describe how to correlate our pause dataset with other datasets e.g. genes, tRNAs and the like.
You should have already calculated pause locations and pause sensitivities by now, or have obtained
these files from the dataset associated with our study or elsewhere (consult `main_workflow.md` if you do not have these files).

The main script you run is called `pipeline_sgm_pauses_winnow_and_analyze.sh`.
Briefly, you need two json files that look like the following excerpts.
The inputs below are all example inputs; you must fill them with locations to real files or realistic parameters.
You can see the comments in that file to see what inputs it needs in detail.

```json
{
  "mod_bam": "/some/dir/mod.sorted.bam",
  "forksense_dir": "/ei/projects/forksense/experiment_XY/",
  "dataset": "20210315_EXP_ONT_XY_abc123",
  "date": "05feb23",
  "output_dir": "/some/dir/pauses_12nov22/analysis_${date}",
  "input_pause_file": "/some/dir/pauses_file",
  "seq_summary_file": "/some/dir/seq_summ.txt",
  "fasta_file": "/some/dir/some_reference.fna",
  "AT_ratio_optional": 0.47
}
```

```json
{
  "pipeline_go_calculate_sensitivity": 1,
  "pipeline_go_calculate_enrichment": 1,
  "analysis_label": "${date}",
  "date": "11mar25",
  "log_file": "/some/log/file/location/log_analysis_${date}.txt",
   "feature_files": [ "/some/feature/file.json",
                      "/some/other/feature/file.json" ]
}
```


The features we are going to correlate with the pauses are in the feature files.
An example is shown below.

```json
{
    "align_by": "head-to-tail",
    "feature": "tRNA",
    "input_bed_file": "/tRNA/locations/zero/size.bed",
    "delete_new_mod_bam_made_on_the_fly": true,
    "genome_size_bp_optional": 10000000,
    "split_by_length_bp": 200,
    "no_jobs_for_all": true,
    "num_split_by_value": 1,
    "split_algorithm_version": "v2",
    "grow_region_bp_if_version_is_v2": 100,
    "fai_file": "/some/dir/some_reference.fna.fai"
}
```