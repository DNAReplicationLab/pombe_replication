We mention two quick scripts to make histograms and scatter plots.

# Make a histogram

We want to make a histogram from an input file `$data_file` that contains binned data in the following format.
We can include other columns and shuffle their order freely.
```text
bin_lo bin_hi value
0 100 0.1
100 200 0.2
```

The columns can be space or tab separated.
As usual, replace the various command-line inputs below suitably.
Some inputs may be optional.

```shell
< $data_file Rscript plot_histogram.R value $plot_png_file $x_axis_label $y_axis_label
```

# Make a scatter plot

We want to make a scatter plot from an input file `$data_file` that contains data in the following format.
We can include other columns and shuffle their order freely.
```text
column_name_x column_name_y colour_column
1 2 a
2 3 a
3 4 a
4 5 b
5 6 b
6 7 b
```

The columns can be space or tab separated.
As usual, replace the various command-line inputs below suitably.
Some inputs may be optional.

```bash
< $data_file Rscript plot_scatter.R column_name_x column_name_y $plot_png_file $x_label $y_label $xlim $ylim $alpha\
    $size colour_column
```