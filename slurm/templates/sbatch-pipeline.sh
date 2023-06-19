#!/bin/bash
# Submit a basic NOMADIC script to Raven
# via Slurm
# JHendry, 2023/06/19

map_jid=$(sbatch --parsable {run_dir}/nomadic-map.sh)
remap_jid=$(sbatch --parsable --dependency=aftercorr:$map_jid {run_dir}/nomadic-remap.sh)

qcbam_jid=$(sbatch --parsable --dependency=aftercorr:$remap_jid {run_dir}/nomadic-qcbam.sh)
targets_jid=$(sbatch --parsable --dependency=aftercorr:$remap_jid {run_dir}/nomadic-targets.sh)

