#!/bin/bash
# Submit a basic NOMADIC script to Raven
# via Slurm
# JHendry, 2023/06/19

bar_jid=$(sbatch --parsable {run_dir}/nomadic-barcode.sh)
expt_jid=$(sbatch --parsable --dependency=afterany:$bar_jid {run_dir}/nomadic-experiment.sh)
