import subprocess
from nomadic.lib.generic import produce_dir


# ================================================================
# ONT barcoding kit availability
#
# ================================================================


BARCODING_KIT_MAPPING = {
    "native24": "EXP-NDB104 EXP-NDB114",
    "native96": "EXP-NBD196",
    "rapid": "SQK-RBK004",
    "pcr": "SQK-PBK004"
}


# ================================================================
# Python interface for guppy_barcoder
#
# ================================================================


def run_guppy_barcode(fastq_input_dir, barcode_kits, output_dir, both_ends, load_module=True):
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
        load_module: bool
            Load guppy usign `module` before executing command. This
            is used when run on BMRC cluster.
    returns
        _ : None
    
    """

    cmd = ""
    if load_module:
        cmd += "module load ont-guppy/5.0.11_linux64 && "

    # Construct command
    cmd += "guppy_barcoder"
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