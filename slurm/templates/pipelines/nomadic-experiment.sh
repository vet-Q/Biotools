#!/bin/bash -l
#SBATCH --job-name=nexpt
#SBATCH --output=logs/nexpt/nexpt-jid%A-%a.out
#SBATCH --error=logs/nexpt/nexpt-jid%A-%a.err
#SBATCH --chdir=./
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=23:00:00


# SETTINGS
expt_dir={expt_dir}
config={config}
bed_path=configs/beds/nomads8.amplicons.bed

# LOAD ANACONDA
module load anaconda/3/.2023.03

# ACTIVATE ENVIRONMENT
conda activate nomadic2-fast

# RUN NOMADIC
nomadic qcbams -e $expt_dir -c $config --overview
nomadic targets -e $expt_dir -c $config --overview
nomadic quickcall -e $expt_dir -c $config -m bcftools -r $bed_path --overview
nomadic quickcall -e $expt_dir -c $config -m clair3_singularity -r $bed_path --overview
nomadic bedcov -e $expt_dir -c $config -r $bed_path --overview
nomadic checkcontam -e $expt_dir -c $config
#nomadic find -e $expt_dir -c $config -m bcftools
#nomadic find -e $expt_dir -c $config -m clair3sing