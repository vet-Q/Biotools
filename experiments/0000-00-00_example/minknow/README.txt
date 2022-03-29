EXAMPLE DATA
2021/08/05, JHendry

This directory includes example fastq files generated
as part of a nanopore run on P. falciparum amplicon data,
sequenced using the native barcoding protocol. I've included
three fastq files from the first three barcodes to allow
easy testing of the pipeline.

Note that the pipeline assumes the fastq data is arranged as follows:

experiments/<expt_name>/<sub_dir>/fastq_pass

The pipeline then produces outputs in:

experiments/<expt_name>/nomadic/<sub_dir>

The motivation for having the subdirectory <sub_dir> is that it 
provides flexibility with respect to rebasecalling or
demultiplexing data from a single experiment. For example,
if there is a guppy update and you want to rebasecall the fast5 files
from an experiment, you can deposit the new fastq files in a new
subdirectory <guppy_vnew>; and then simple run nomadic on the 
subdirectory.

