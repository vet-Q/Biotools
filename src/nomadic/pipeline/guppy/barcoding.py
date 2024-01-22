import subprocess
from nomadic.lib.generic import produce_dir


# ================================================================
# ONT barcoding kit availability
#
# ================================================================


BARCODING_KIT_MAPPING = {
    "native24": '"EXP-NBD104 EXP-NBD114"',
    "native96": "EXP-NBD196",
    "rapid": "SQK-RBK004",
    "pcr": "SQK-PBK004",
    "R10_native96": "SQK-NBD114-96",
    "R10_rapid96": "SQK-RBK114-96"
}


# ================================================================
# Python interface for guppy_barcoder
#
# ================================================================


def run_guppy_barcode(
    fastq_input_dir: str,
    barcode_kits: str,
    output_dir: str,
    both_ends: bool = False,
    strict: bool = False,
    recursive: bool = False,
    dry_run: bool = False,
) -> None:
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

    # Construct command
    cmd = "guppy_barcoder"
    cmd += " --device 'cuda:0'"
    cmd += f" --barcode_kits {barcode_kits}"
    cmd += f" --input_path {fastq_input_dir}"
    if recursive:
        cmd += " ---recursive"
    if both_ends:
        cmd += " --require_barcodes_both_ends"
    if strict:
        cmd += " --min_score_barcode_front 90"
        cmd += " --min_score_barcode_rear 90"
    cmd += f" --save_path {output_dir}"
    cmd += " --enable_trim_barcodes"
    cmd += " --compress_fastq"
    cmd += " --disable_pings"

    if dry_run:
        print(cmd)
        return

    # Run
    subprocess.run(cmd, shell=True, check=True)
