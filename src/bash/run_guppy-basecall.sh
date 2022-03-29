#/bin/bash

# LOOAD GUPPY
module load ont-guppy/5.0.11_linux64

# PARSE CLI
echo "Parsing user inputs..."
cli=$@
while [[ $# -ge 1 ]]
do
  key=$1
  case $key in
  	-e)
		expt_dir=$2
		shift
		;;
	-m)
		method=$2
		shift
		;;
	*)
		;;
  esac
  shift
done

# Make output directory
input_dir=$(echo $expt_dir/minknow)
output_dir=$(echo $expt_dir/guppy/$method)
mkdir -p $output_dir

# Define configuration
config_fn=$(echo "dna_r9.4.1_450bps_"$method".cfg")

# Print to stdout
echo "Experiment dir:" $expt_dir
echo "Input dir:" $input_dir
echo "Basecalling method:" $method
echo "  Configuration:" $config_fn
echo "Output dir:" $output_dir

# Run guppy
echo "Running guppy basecaller.."
guppy_basecaller \
--device "cuda:0" \
--min_qscore 8 \
--compress_fastq  \
--config $config_fn \
--input_path $input_dir \
--recursive \
--save_path $output_dir \
--disable_pings
echo "Done."
echo ""
