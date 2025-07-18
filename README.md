# pombe_replication
---
## About the repository

This repository contains the scripts used for data analysis and plotting in the fission yeast replication dynamics manuscript "Single-molecule sequencing maps replication dynamics across the fission yeast genome, including centromeres". 

mod.bam and text files (i.e., initiation sites, termination sites, pause sites, leftward fork ratio, reference genome sequences, annotations, replication timing data) are available from the Zenodo repository under the DOI [10.5281/zenodo.15365203].


- Instructions, scripts and model used to detect BrdU in nanopore reads, and to call replication dynamics (initiation, termination, fork directionalities and pauses) are in the repository [fork_arrest](https://github.com/DNAReplicationLab/fork_arrest).

- Scripts used to produce most of the figures of the cited manuscript are in the folder `pombe_figures_scripts/`. 

## Structure of the repository

- `pombe_replication/`
  - `fork_arrest_scripts/`  
    → DNAscent pipeline and pause detection pipeline from repository [fork_arrest](https://github.com/DNAReplicationLab/fork_arrest)
  - `pombe_figures_scripts/`  
    → Scripts used in the fission yeast replication dynamics manuscript
  - `README.md`


## Instructions
Please follow these instructions in order to run the code.

1. Download the `fork_arrest_scripts/` repository.
2. Download the scripts from `pombe_figures_scripts/` files and save them in the `scripts/` folder from the `fork_arrest_scripts/` repository.
3. Within `scripts/`, create the `config.sh`, `config.py`, and `config.R` files as described in the file `DNAReplicationLab/fork_arrest/README.md`.
4. Download the files from the Zenodo repository under the DOI [10.5281/zenodo.15365203] and save them in `input_dir/`.
5. Index the `.mod.bam` files using the following command:  
   ```bash
   samtools index -M input_dir/*.mod.bam
6. Create a directory where the output files will be stored: output_dir/.
7. Run the code from `scripts/`.
For example:
```bash
   cd fork_arrest/scripts/
   bash Fig1B-D_rainplots.sh /path/to/input_dir /path/to/output_dir
