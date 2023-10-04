#!/bin/bash -l
#SBATCH --job-name=dorado
#SBATCH --output=logs/dorado-jid%A-%a.out
#SBATCH --error=logs/dorado-jid%A-%a.err
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

echo "Run dorado on MPCDF."
echo "Currently in directory:"`pwd`

# Settings
EXPT_DIR={expt_dir}
CONFIG={config}  # here, not used
ACCURACY={basecalling_method}

# NB: you have to download the models first
BASECALL_MODEL="dorado_models/dna_r10.4.1_e8.2_400bps_"$ACCURACY"@v4.2.0"

POD5_DIR=$EXPT_DIR"/minknow"
FASTQ_DIR=$EXPT_DIR"/dorado/"$ACCURACY
FASTQ_PATH=$FASTQ_DIR/"dorado.called.fastq"

mkdir -p $OUTPUT_DIR
mkdir -p logs

echo "-------------------------------------------------------------------"
echo "Preparing to run Dorado"
echo "  Experiment Dir.: "$EXPT_DIR
echo "  Accuracy: "$ACCURACY
echo "  Basecall model: "$BASECALL_MODEL
echo "  POD5 dir.: "$POD5_DIR
echo "  Output dir.: "$FASTQ_DIR
echo "  Output FASTQ: "$FASTQ_PATH

echo "Runninng..."
dorado basecaller $BASECALL_MODEL $POD5_DIR \
--device 'cuda:0' \
--min-qscore 10 \
--emit-fastq \
--recursive > $FASTQ_PATH
echo "Done."
echo ""

echo "-------------------------------------------------------------------"
echo "Preparing to demultiplex with Guppy"
echo "Currently in directory:"`pwd`

DEMUX_DIR=$FASTQ_DIR"/single_end_strict"
BARCODE_KIT={barcoding_strategy}

mkdir -p $DEMUX_DIR

echo "Preparing to run guppy"
echo "  FASTQ dir.: "$FASTQ_DIR
echo "  Demultiplex dir.: "$DEMUX_DIR
echo "  Barcode kit: "$BARCODE_KIT

echo "Runninng..."
guppy_barcoder \
--device 'cuda:0' \
--barcode_kits $BARCODE_KIT \
--recursive \
--compress_fastq \
--enable_trim_barcodes \
--trim_adapters \
--input_path $FASTQ_DIR \
--save_path $DEMUX_DIR \
--min_score_barcode_front 90 \
--min_score_barcode_rear 90 \
--disable_pings
echo "Done."
echo ""

echo "-------------------------------------------------------------------"