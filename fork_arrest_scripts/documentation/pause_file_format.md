The pause file is a tab-separated value format. Each row corresponds to one single-molecule
replisome pause detection.

1. **detectIndex**: A string in the format `readID_contig_start_end_orn_dirn_startFork_endFork`.
   - **readID**: Identifier for the read.
   - **contig**: Contig name, which can contain underscores.
   - **start**: Start position of the read on the reference genome (integer).
   - **end**: End position of the read on the reference genome (integer).
   - **orn**: Orientation, either `fwd` (forward) or `rev` (reverse).
   - **dirn**: Direction, either `L` (left) or `R` (right).
   - **startFork**: Start position of the fork on the reference genome (integer).
   - **endFork**: End position of the fork on the reference genome (integer).

2. **pauseSite**: The position along the reference genome where the pause occurs.

3. **keep***: Optional boolean columns used to filter the data. If any of the `keep*` columns is `False`, the pause is discarded
(if the `--discardUnkeptPauses` option is used in our various scripts e.g. `convert_pause_file_to_bed.py`).

The pause file can be convert to the bed format using `convert_pause_file_to_bed.py`.
Please have a look at the various options there, for e.g.: you can use the fork direction to set the strand in the bed file.
Convertion to a bed file is useful for downstream analysis.