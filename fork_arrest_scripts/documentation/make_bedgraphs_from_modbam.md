We discuss here how to make bedgraphs of data from modBAM files.

# Make bedgraphs of raw or windowed data from one read

```shell
# raw data
python get_raw_data_from_modBAM.py $mod_bam $read_id chrVI 50000 60000 |\
	  awk '{print $1 " chrVI " $2 " " $3}' |\
	  sed '1i\index contig start val' |\
	  python make_bedgraphs.py --dir $plot_dir --ncol 3
```

To window data, one needs more parameters like window size (in thymidines) and BrdU thresholds,
which have been set to 300 and 0.5 respectively in the snippet below. Hence, each window has 300
thymidines, each of which is considered as modified or unmodified depending on whether the probability
is above or below 0.5.

```shell
# windowed data
python get_raw_data_from_modBAM.py $mod_bam $read_id chrVI 50000 60000 |\
	  sed '1i\detectIndex posOnRef val' |\
	  python get_mean_brdU_window.py --window 300 --thres 0.5 --halfOpen |\
	  awk '{print $0 " chrVI"}' |\
	  sed '1i\index val start end contig' |\
	  python make_bedgraphs.py --dir $plot_dir --ncol 4 --suffix _win
```

# Make bedgraphs of raw or windowed data from many reads

Please read the one read versions above first.

Prepare a regions file with five space-separated columns with the headers
`read_id contig start end alt_read_id`. The columns can be in any order.
Ensure that alt read id contains the read id and the contig. An example
for alt_read_id is `readid_contig_start_end`. Then, the parameter `--implicitContigFormat`
can be set to `_2` which means that contig information is the second bit of information
in the index using the delimiter `_`.

```shell
# raw data
< $regions_file python get_raw_data_from_modBAM.py --piped-regions --alt-read-id-column $modbam |\
	  sed '1i\index start val' |\
	  python make_bedgraphs.py --dir $plot_dir --ncol 3  --implicitContigFormat _2
```

```shell
# windowed data
< $regions_file python get_raw_data_from_modBAM.py --piped-regions --alt-read-id-column $modbam |\
	  sed '1i\detectIndex posOnRef val' |\
	  python get_mean_brdU_window.py --window 300 --thres 0.5 --halfOpen |\
	  sed '1i\index val start end' |\
	  python make_bedgraphs.py --dir $plot_dir --ncol 4  --implicitContigFormat _2 --suffix _win
```