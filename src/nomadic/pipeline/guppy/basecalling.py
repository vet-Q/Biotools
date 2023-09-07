import subprocess


# ================================================================
# Define available basecalling methods
#
# ================================================================


FLOW_CELLS = {"R9": "dna_r9.4.1_450bps", "R10": "dna_r10.4.1_e8.2_400bps"}
BASECALL_METHODS = ["sup", "hac", "fast"]


# ================================================================
# Python interface for guppy_basecaller
#
# ================================================================


def run_guppy_basecaller(
    input_dir: str,
    output_dir: str,
    flow_cell: str = "R10",
    method: str = "sup",
    min_qscore: int = 12,
    dry_run: bool = False,
) -> None:
    """
    Run guppy basescaller

    Notes
    -----
    1. We are assuming `guppy_basecaller` is available from the
    command line

    2. Underneath the `output_dir`, `guppy` will write...
        /pass
        /fail
    ...on the basis of the `min_qscore`. For R10 data, a minimum
    Q-score of 12 probably makes sense, it gets rid of reads that might
    otherwise be misclassified or cause issues in downstream analyses.

    3. The default input directory is just `minknow` and `--recursive` is used.
    This will work regardless of whether or not the pod5 are in `pod5` or
    `pod5_pass` and `pod5_fail`.

    """

    # Configuration path depends on basecallig method
    if flow_cell not in FLOW_CELLS:
        raise ValueError(
            f"`flow_cell` {flow_cell} must one of {', '.join(FLOW_CELLS)}."
        )

    if method not in BASECALL_METHODS:
        raise ValueError(
            f"basecalling `method` {method} must be one of {', '.join(BASECALL_METHODS)}."
        )

    config_path = f"{FLOW_CELLS[flow_cell]}_{method}.cfg"

    # # Optionally load module (for BMRC)
    # cmd = ""
    # if load_module:
    #     cmd += "module load ont-guppy/5.0.11_linux64 && "

    # Construct command
    cmd = "guppy_basecaller"
    cmd += " --device 'cuda:0'"
    cmd += f" --min_qscore {min_qscore}"
    cmd += f" --config {config_path}"
    cmd += f" --input_path {input_dir}"
    cmd += f" --save_path {output_dir}"
    cmd += " --compress_fastq"
    cmd += " --recursive"
    cmd += " --disable_pings"

    if dry_run:
        print(cmd)
        return

    # Run
    subprocess.run(cmd, shell=True, check=True)
