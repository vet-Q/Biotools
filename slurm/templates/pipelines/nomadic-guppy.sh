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

# SETTINGS
# Directories
POD5_DIR=$expt_dir"/minknow/pod5_pass"
FASTQ_DIR=$expt_dir"/guppy/hac/fastq/pass"
DEMUX_FASTQ_DIR=$expt_dir"/guppy/hac/single_end"

mkdir -p $FASTQ_DIR
mkdir -p $DEMUX_FASTQ_DIR

# Parameters
BASECALL_SETTINGS="dna_r10.4.1_e8.2_400bps_hac.cfg"
BARCODE_KIT="SQK-NBD114-96"

# LOAD ANACONDA & ACTIVATE
module load anaconda/3/.2023.03
conda activate nomadic2-fast

# RUN BASECALLING
echo "Runninng basecalling..."
guppy_basecaller \
--device 'cuda:0' \
--min_qscore 8 \
--compress_fastq \
--config $BASECALL_SETTINGS \
--input_path $POD5_DIR \
--recursive \
--save_path $FASTQ_DIR \
--disable_pings
echo "Done."
echo ""

# RUN BARCODING
echo "Runninng barcoding..."
guppy_barcoder \
--device 'cuda:0' \
--recursive \
--compress_fastq \
--enable_trim_barcodes \
--trim_adapters \
--input_path $FASTQ_DIR \
--save_path $DEMUX_FASTQ_DIR \
--disable_pings \
#--require_barcodes_both_ends
echo "Done."
echo ""

