#!/bin/bash -l
#SBATCH --job-name=nguppy
#SBATCH --output=logs/nguppy/nguppy-jid%A-%a.out
#SBATCH --error=logs/nguppy/nguppy-jid%A-%a.err
#SBATCH --chdir=./
#SBATCH --ntasks=1
# --------------------
# GPU Commands
# --------------------
#SBATCH --constraint="gpu"
#SBATCH --gres=gpu:a100:1
#SBATCH --cpus-per-task=18
#SBATCH --mem=125000
#
#SBATCH --time=24:00:00


# SETTINGS
expt_dir={expt_dir}
config={config}
basecalling_method={basecalling_method}
barcode_kit={barcoding_strategy}

# LOAD ANACONDA
module load anaconda/3/.2023.03

# ACTIVATE ENVIRONMENT
conda activate nomadic2-fast

# RUN NOMADIC
nomadic basecall \
-e $expt_dir \
-f R10 \
-m $basecalling_method \
-q 12

nomadic barcode \
-e $expt_dir \
-m $basecalling_method \
-k $barcode_kit {barcode_strict} {require_both_ends}

