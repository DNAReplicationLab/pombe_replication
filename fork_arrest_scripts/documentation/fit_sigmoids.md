After getting analogue densities in mod BAM files and replication track annotations
from forkSense, you can now fit sigmoids per fork (consult `main_workflow.md` if you do not have these files).
Use the following command:

```bash
mod_bam=/path/to/mod.bam
forkSense_directory=/path/to/forkSenseOverallBedgraphs
op_dir=/path/to/output/forkLen10AlignLen30/modelFitResults/fork_fit_1000
n_forks=1000

dateStr=09apr25

mkdir -p "$op_dir"

for tag in {A..Z}; do
  slurm_out="$op_dir"/output_"$dateStr"_"$tag".out
  slurm_err="$op_dir"/output_"$dateStr"_"$tag".err
  sbatch -p ei-short --time=11:59:59 -o "$slurm_out" -e "$slurm_err" \
    subsample_forks_fit_cg_model_and_plot.sh "$mod_bam" "$forkSense_directory" "$op_dir" "$dateStr"_"$tag" "$n_forks"
done
```

This will produce ~26 fits in that folder (I say approximately because there is a chance some jobs may fail).
You can then average best-fit parameters across all runs to get the parameters of the reference sigmoid.
There are various assumptions in the model fitting; please see the script to know what they are.
One in particular is there are many more than 26000 'long' forks in our dataset as we sample 1000 forks independently
at a time and do 26 runs with the hope that the same fork is not picked repeatedly.
Also we scan over a grid of candidate values for the low, high, and width parameters of the sigmoid.
If your experiment is such that you get best-fits at the edge of or beyond the grid, then you've
to change the bounds of the grid manually in the script.