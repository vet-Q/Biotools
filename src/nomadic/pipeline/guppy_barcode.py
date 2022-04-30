import subprocess
from nomadic.lib.generic import produce_dir

# ================================================================
# Parameters
#
# ================================================================


ONLY_PASS = True  # only demultiplex .fastq that pass quality control
BARCODING_KIT_MAPPING = {
    "native24": "EXP-NDB104 EXP-NDB114",
    "native96": "EXP-NBD196",
    "rapid": "SQK-RBK004",
    "pcr": "SQK-PBK004"
}


# ================================================================
# Parameters
#
# ================================================================


def run_guppy_barcode(fastq_input_dir, barcode_kits, output_dir, both_ends):
    """
    Run guppy demultiplexing

    params:
        fastq_input_dir: str
            Path to directory containing .fastq files to demultiplex. Using
            `--recursive` so can be across multiple sub-directories.
        barcode_kits: str
            A space-separated string of valid ONT barcoding kits.
        output_dir: str
            Output directory; must already exist.
        both_ends: bool
            Should demultiplexing require barcodes on *both* ends?
    returns
        _ : None
    
    """

    # Construct command
    cmd = "guppy_barcoder"
    cmd += " --device 'cuda:0'"
    cmd += " --compress_fastq"
    cmd += " --trim_barcodes"
    cmd += f" --barcode_kits {barcode_kits}"
    cmd += f" --input_path {fastq_input_dir}"
    cmd += " --recursive"
    cmd += f" --save_path {output_dir}"
    cmd += " --disable_pings"
    if both_ends:
        cmd += " --require_barcodes_both_ends"
    
    # Run
    subprocess.run(cmd, shell=True, check=True)


# ================================================================
# Main script, run from `cli.py`
#
# ================================================================


def main(expt_dir, basecalling_method, barcoding_strategy, both_ends):
    """
    Run guppy demultiplexing on .fastq files

    """
    # LOAD GUPPY
    subprocess.run("module load ont-guppy/5.0.11_linux64", shell=True, check=True)

    # SELECT KIT
    barcode_kits = BARCODING_KIT_MAPPING[barcoding_strategy]

    # CREATE OUTPUT DIRECTORY
    input_dir = f"{expt_dir}/guppy/{basecalling_method}"
    if ONLY_PASS:
        fastq_input_dir = f"{input_dir}/pass"
    output_dir = produce_dir(input_dir, 'both_ends' if both_ends else 'single_end')

    # PRINT TO STDOUT
    print("Inputs")
    print(f"  Experiment dir.: {expt_dir}")
    print(f"  Basecalling method: {basecalling_method}")
    print(f"  Input directory: {input_dir}")
    if both_ends:
        print(f"  Requiring both ends to be barcoded.")
    else:
        print("  Requiring only one end to be barcoded.")
    print(f"  Barcode kits: {barcode_kits}")
    print(f"  Will run with standard configuration.")
    print(f"  Output directory: {output_dir}")
    print("Done.")
    print("")

    # RUN GUPPY
    print("Running guppy barcoder...")
    run_guppy_barcode(
        fastq_input_dir=fastq_input_dir,
        barcode_kits=barcode_kits,
        output_dir=output_dir,
        both_ends=both_ends
    )
    print("Done.")
    print("")