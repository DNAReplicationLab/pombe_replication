# Instructions to manipulate pod5 files

## Introduction

Historically nanopore raw data (timepoints of electric currents) was stored in fast5 files.
There may be one read per file or multiple reads per file. These files are not human-readable.

Now, there is a faster format called pod5 and fast5 is being phased out.
Thus, we need tools to help us interconvert between fast5 and pod5 files.
Pod5 files are not human-readable as well.

It is important to note that modification calls produced from either of these file formats are stored by us
in the universal mod BAM format. So the format of the raw data is not important for any tool that receives
modification data in the universal mod BAM format as input.

## The pod5 package is used to interact with pod5 files

ONT has developed a package called [pod5](https://github.com/nanoporetech/pod5-file-format) to process
pod5 files and convert fast5 files to pod5. The package has already been installed in the Earlham HPC and can be loaded
in an interactive job or in a script using the command `source load_package.sh -pod5` executed from the root folder of
this repository.

### Caution about commands below

* Although we have tested most of the commands below, we cannot guarantee that they will work in all cases.
* If a command does not work, please let us know.
* If a command is in conflict with the official documentation, then the official documentation should be followed.
* It is up to the user to allocate enough resources (memory, runtime etc.) to the command. It is a good rule of
  thumb that processing of large files should be done on an HPC using for example a slurm job with the appropriate
  resources.

### Example pod5 commands to convert fast5 to pod5 and vice versa

The following commands convert a fast5 file to a pod5 file and vice versa.

```bash
pod5 convert from_fast5 $fast5_file -o $pod5_file
pod5 convert to_fast5 $pod5_file -o $fast5_folder # NOTE that output is sent to a folder
```

* The commands `pod5 convert from_fast5 -h`, `pod5 convert to_fast5 -h` can be used to
  get more information about the commands.
* Use the flag `--recursive` and a directory instead of an individual file as input to recursively convert all files
  in that directory. This command may even descend into subdirectories; we do not know.
* Use the parameter `-t $num_threads` to specify the number of threads to use. The default is 4.
* By default, only 4000 reads are stored per fast5 file. So if you convert one pod5 file with many reads, you may
  get many fast5 files. If you want to change this behaviour, please use the flag `--file-read-count $num_reads`.
* If you want to overwrite a pre-existing file(s), please use the flag `-f`.

### Example pod5 command to retain specific reads

This command reads all the pod5 files in the folder `"$pod5_folder"` and retains only the reads that are
listed in the file `"$some_folder"/read_ids.txt` (text file should contain one read per line and no header or
column name). The output is stored in the file `"$output_folder"/subset.pod5`.

```bash
pod5 filter "$pod5_folder"/* -i "$some_folder"/read_ids.txt \
    -o "$output_folder"/subset.pod5
```

### Other pod5 commands of interest

* `pod5 merge` can be used to merge multiple pod5 files into one pod5 file.
* `pod5 inspect`/`pod5 view` can be used to extract summary information from a pod5 file.