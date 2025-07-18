We discuss how to make a genomic annotation image with basic shapes like arrows, circles etc.
that correspond to different features.

# Make an image of genomic annotations, optionally to accompany plots of data vs genomic coordinates

Example plot of fake data below, made from `sample_files/sample_genetic_elements.json`
![Sample genomic elements plot](../sample_files/sample_genetic_elements.svg "Sample genomic elements plot")

Command to produce the plot

```shell
python plotting_and_short_analyses/make_genomic_annotation_image.py\
  --annotate sample_files/sample_genetic_elements.json --op sample_files/sample_genetic_elements.svg
```

A small snippet of the json file is shown below for illustrative purposes. The image is of width `1000`,
where x coordinates of `140` and `700` map to genomic coordinates of `120` and `165`. A reference line
runs from `120` to `160`, and contains a gene from `120` to `130` and an active origin at `120`. Objects
are scaled vertically by a length corresponding to `1 kb` on the x-axis through the parameter `scale_factor`.
Many more symbols and options are illustrated in the json file.

```json
{
  "scale_factor": 1,
  "img_width": 1000,
  "mapping_coordinates_1": [
    120,
    140
  ],
  "mapping_coordinates_2": [
    165,
    700
  ],
  "genomic_elements": [
    {
      "type": "reference_line",
      "start": 120,
      "stop": 160
    },
    {
      "type": "gene",
      "start": 120,
      "stop": 130,
      "modifier": "emphasize"
    },
    {
      "type": "origin",
      "position": 120,
      "modifier": "active"
    }
  ]
}
```

A similar visualization can also be produced using our fork of the DNA features viewer package.
The command and visualization are below; the figure is made from `sample_files/sample_genetic_elements_dfv.json`
![Sample genomic elements plot](../sample_files/sample_genetic_elements_dfv.svg "Sample genomic elements plot from DFV")

Command to produce the plot

```shell
python plotting_and_short_analyses/make_genomic_annotation_image.py\
  --annotate sample_files/sample_genetic_elements_dfv.json --op sample_files/sample_genetic_elements_dfv.svg --use-dna-features-viewer
```