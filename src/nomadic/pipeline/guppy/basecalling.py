import subprocess


# ================================================================
# Define available basecalling methods
#
# ================================================================


BASECALL_METHODS = ["hac", "fast"]


# ================================================================
# Python interface for guppy_basecaller
#
# ================================================================


def run_guppy_basecaller(input_dir, output_dir, method="hac", load_module=True):
    """
    Run guppy basescaller

    """

    # Configuration path depends on basecallig method
    if method not in BASECALL_METHODS:
        raise ValueError(
            f"`basecall_method` {method} must be one of {', '.join(BASECALL_METHODS)}."
        )
    config_path = f"dna_r9.4.1_450bps_{method}.cfg"

    # Optionally load module (for BMRC)
    cmd = ""
    if load_module:
        cmd += "module load ont-guppy/5.0.11_linux64 && "

    # Construct command
    cmd += "guppy_basecaller"
    cmd += " --device 'cuda:0'"
    cmd += " --min_qscore 8"
    cmd += " --compress_fastq"
    cmd += f" --config {config_path}"
    cmd += f" --input_path {input_dir}"
    cmd += " --recursive"
    cmd += f" --save_path {output_dir}"
    cmd += " --disable_pings"

    # Run
    subprocess.run(cmd, shell=True, check=True)
