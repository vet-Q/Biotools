#!/bin/bash -l
#SBATCH --job-name=ncalldown
#SBATCH --output=logs/ncalldown/ncalldown-jid%A-%a.out
#SBATCH --error=logs/ncalldown/ncalldown-jid%A-%a.err
#SBATCH --chdir=./
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array {array_str}
#SBATCH --mem=16GB
#SBATCH --time=23:00:00


# SETTINGS
expt_dir={expt_dir}
config={config}

# LOAD ANACONDA
module load anaconda/3/.2023.03

# ACTIVATE ENVIRONMENT
conda activate nomadic2-fast

# To run Clair3 via Singularity, need to load singularity module
module load singularity/link2apptainer

# Run Clair3 with downsampling
nomadic call \
-e $expt_dir \
-c $config \
-m clair3sing \
--downsample \
-r 10 -r 20 -r 30 -r 40 -r 50 -r 60 -r 70 -r 80 -r 90 -r 100 \
-i 10 \
-b $SLURM_ARRAY_TASK_ID