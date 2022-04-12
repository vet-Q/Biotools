NOMADIC EXAMPLE DATA
2022/04/12, J. Hendry


### Basic Setup

`nomadic` assumes that sequencing data and metadata are
arranged in a particular fashion. Follow the steps below
to prepare your data correctly.

First, give your experiment a name <expt_name>. Create
a directory:

mkdir experiments/<expt_dir>

Create two sub-directories:

mkdir experiments/<expt_dir>/minknow
mkdir experiments/<expt_dir>/metadata

Inside of `experiments/<expt_dir>/metadata`, you should
place a `.csv` file containing metadata about your samples.
This file must contain the following columns:

barcode 
    A list of barcodes used during library prep.
    e.g., barcode01, barcode02, &c.

sample_id
    A unique identifier for each sample.
    e.g. samp001, samp002, &c.

Other columns can be included.  See:

experiments/0000-00-00_example/metadata/sample_info.csv 

...for an example.

Inside of `experiments/<expt_dir>/minknow`, you should
deposit the output folder from MinKNOW. This should include
a `/fastq_pass` directory containing .fastq files.

Note that this arrangement should match the `[Experiment]`
section of your configuration file, e.g.

```
[Experiment]
metadata = metadata/sample_info.csv
basecalling = minknow
```


### Alternative sub-directories

The same set of `.fast5` files can be basecalled multiple times,
with differet basecallers, versions of `guppy`, &c. `nomadic`
supports multiple basecalling runs on the same experiment. Simply
create a new sub-directory, e.g.

mkdir experiments/<expt_dir>/new_basecalling

...and modify the your configuration file in /configs to indicate
that .fastq files are stored here. 

```
[Experiment]
metadata = metadata/sample_info.csv
basecalling = new_basecalling
```

