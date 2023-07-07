#!/bin/bash
# Submit a basic NOMADIC script to Raven
# via Slurm
# JHendry, 2023/06/19

gup_jid=$(sbatch --parsable {run_dir}/nomadic-guppy.sh)
bar_jid=$(sbatch --parsable --dependency=afterany:$gup_jid {run_dir}/nomadic-barcode.sh)
expt_jid=$(sbatch --parsable --dependency=afterany:$bar_jid {run_dir}/nomadic-experiment.sh)
