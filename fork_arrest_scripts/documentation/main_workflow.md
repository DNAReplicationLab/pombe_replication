There are three main pipelines in this repository:
- Consult `fast5_to_dnascent_pipeline.md` to learn how to go from nanopore currents to
  analogue calls and replication-track annotations.
- Consult `fit_sigmoids.md` to fit sigmoids to the identified tracks to get a reference
  sigmoid.
- Consult `pause_finding_pipeline.md` to use the reference sigmoid to obtain pauses.
- Consult `pause_analysis_pipeline.md` to correlate pause sites with features on the genome.
- Consult `pause_report_in_region.md` to correlate pause sites with one bed file.

The above steps can be used if you want to perform pause calling anew.
If you already have some or all of these files e.g. from a previous run or from the dataset associated
with this study, then you can avoid some or all of these steps.

We use a simpler pipeline to find pauses in the rDNA.
The details of this are in `pause_finding_pipeline.md`.
To calculate pauses in the rDNA, you do not need to perform sigmoid fitting i.e. you do not need to
consult `fit_sigmoids.md`.