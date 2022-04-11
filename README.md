<p align="center"><img src="images/nomadic_logo-01.png" width="500"></p>

## Overview

NOMADIC is a collection of python scripts for analysing long read data generated from *Plasmodium falciparum* malaria. The scripts are assembled into a pipeline that can be run locally or on an Oracle Grid Engine cluster. 

At present, the focus is nanopore data generated for a panel of target amplicons.

## Install
First, we build a virtual environment including all dependencies using `conda`:

```
conda update conda
conda env create
conda activate nomadic2
```

Next, we locally install the python package `nomadic2` into this virtual environment using `pip`:

```
pip install .
```

## Quick start
First, output data from `minknow` and a metadata file must be arranged into the `experiments` directory. See `experiments/0000-00-00_example` for an example. Then, the full pipeline can be run by invoking:

```
nomadic runall -e experiments/0000-00-00_example -c configs/default.ini
```

Where `0000-00-00_example` would be replaced with your experiment of interest. The configuration file contains information about basecalling settings, target amplicons and mutations. It can be changed to suit your assay.

## Usage

A full list of commands can be found by typing `nomadic`:

```
Usage: nomadic [OPTIONS] COMMAND [ARGS]...

  NOMADIC: A pipeline for analysis of malaria long-read data

Options:
  --help  Show this message and exit.

Commands:
  bmrc     Build BMRC pipeline submission.
  map      Map to P.f. reference.
  qcbams   QC analysis of .bam files.
  remap    Map unmapped reads to H.s.
  runall   Run complete pipeline.
  targets  Analyse amplicon targets.
  
```

