#/bin/bash

# LOAD GUPPY
module load ont-guppy/5.0.11_linux64

# PARSE CLI
echo "Parsing user inputs..."
cli=$@
while [[ $# -ge 1 ]]; do
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
	-k)
		kit=$2
		shift
		;;
	-b)
		both_ends=T
		shift
		shift
		;;
	*)
		;;
  esac
  shift
done

# Prepare barcoding kit
if [[ $kit = native ]]; then
  # EXP-NDB104 ; for 1-12
  # EXP-NDB114 ; for 13-24
  # EXP-NBD196 ; for 1-96
  kit_args="EXP-NBD104 EXP-NBD114 SQK-RBK004"
  # kit_args="EXP-NBD196"
elif [[ $kit = pcr ]]; then
  kit_args="SQK-PBK004"
elif [[ $kit = rapid ]]; then
  # SQK-RBK004 ; for rapid barcodes
  kit_args="SQK-RBK004"
else
  echo "Invalid kit selection."
  echo "Choose from {native, pcr, rapid}."
  exit 1
fi

# Make output directory
input_dir=$(echo $expt_dir/guppy/$method)
only_pass=T
if [[ $only_pass = T ]]; then
  fastq_input_dir=$(echo $input_dir/pass)
fi
if [[ $both_ends = T ]]; then
  output_dir=$(echo $input_dir/both_ends)
else
  output_dir=$(echo $input_dir/single_end)
fi
mkdir -p $output_dir

# Print to stdout
echo "Experiment dir:" $expt_dir
echo "Basecalling method:" $method
echo "Input dir:" $input_dir
if [[ $both_ends = 'T' ]]; then
  echo "Requiring both ends to be barcoded."
else
  echo "Requiring one end to be barcoded."
fi
echo "Barcode kit(s):" $kit_args
echo "Will run with standard configuration."
echo "Output dir:" $output_dir

# Run guppy
echo "Running barcoding..."
guppy_barcoder \
--device "cuda:0" \
--compress_fastq  \
--trim_barcodes \
--barcode_kits "$kit_args" \
--input_path $fastq_input_dir \
--recursive \
--save_path $output_dir \
--disable_pings \
${both_ends:+--require_barcodes_both_ends}
echo "Done."
echo ""

