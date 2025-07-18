It is useful to have a document summarizing statistics measured from a nanopore experiment with modified bases.
We have a tool called `dnascent_qc.sh` that performs just that.
A sample output file `sample_dnascent_qc.sh` is provided in the documentation (please do not assume this is real data).
A sample command is given below, which can also be run in a bash job with sufficient memory.
Substitute variables suitably to point to the correct files and directories.
If you do not have some files like the `$forkSense_dir` for example, you may be able to leave them blank;
please consult the header of the script for how to do this.

```shell
sbatch dnascent_qc.sh $sequencing_summary $mod_bam $forksense_dir $op_dir
# needs sequencing summary file, .mod.bam generated from dnascent detect,
# directory where forksense puts forks, origins, terminations
# output directory where the qc report will be sent with a suitable name like dnascent_qc.pdf
```