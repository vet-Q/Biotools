import os
import click

from pathlib import Path

from nomadic.lib.generic import produce_dir
from nomadic.lib.parsing import build_parameter_dict
from nomadic.pipeline.cli import experiment_options, barcode_option
from nomadic.pipeline.guppy.barcoding import BARCODING_KIT_MAPPING
from nomadic.pipeline.guppy.basecalling import BASECALL_METHODS


GUPPY_PIPELINE = Path("slurm/templates/pipelines/nomadic-guppy.sh")
BAR_PIPELINE = Path("slurm/templates/pipelines/nomadic-barcode.sh")
EXPT_PIPELINE = Path("slurm/templates/pipelines/nomadic-experiment.sh")
DWNSAMP_PIPELINE = Path("slurm/templates/pipelines/nomadic-call-downsample.sh")
SUBMIT_PIPELINE = Path("slurm/templates/nomadic-submit.sh")
RUNS_DIR = Path("slurm/runs")


def load_format_write(input_file: Path, output_file: Path, **kwargs) -> None:
    """Load a file that has named formating fields, e.g. {job_name}, format it, and write"""
    
    input_str = "".join(open(input_file, "r").readlines())
    output_str = input_str.format(**kwargs)
    with open(output_file, "w") as output:
        output.write(output_str)



@click.command(short_help="Generate a NOMADIC submission script for Raven.")
@experiment_options
@barcode_option
@click.option(
    "-m",
    "--basecalling_method",
    type=click.Choice(BASECALL_METHODS),
    default="hac",
    help="Basecalling method; super, high or fast accuracy.",
)
@click.option(
    "-k",
    "--barcoding_strategy",
    type=click.Choice(BARCODING_KIT_MAPPING),
    default="R10_native96",
    help="Barcoding strategy, human-readable names that map to ONT kits.",
)
@click.option(
    "-b",
    "--both_ends",
    is_flag=True,
    default=False,
    help="Require both ends to have barcode?",
)
@click.option(
    "-s",
    "--strict",
    is_flag=True,
    default=False,
    help="Only classify barcodes if alignment score exceeds a strict threshold.",
)
def main(expt_dir: str, config: str, barcode: str, basecalling_method: str, barcoding_strategy: str, both_ends: bool, strict: bool):
    """
    This script will produce a bash script that will launch the NOMADIC pipeline
    for a specific experiment on Raven

    Array jobs are inferred from the metadata file, and dependencies are set
    such that the whole pipeline will be queued when the script is run

    """
    
    # Load parameters
    params = build_parameter_dict(expt_dir, config, barcode)
    
    # Prepare array string
    print(f"Found {len(params['barcodes'])} barcodes to process...")
    array_str = ",".join([str(int(b[-2:])) for b in params["barcodes"]]) # assuming <name><val>

    # Prepare output directory
    expt_dir = Path(expt_dir)
    expt_name = expt_dir.name
    run_count = sum([1 for d in os.listdir(RUNS_DIR) if d.startswith(expt_name)])  # we may create multiple submission /w diff. settings
    output_dir = Path(produce_dir(RUNS_DIR, f"{expt_name}-r{run_count}"))

    # Formatting arguments for guppy
    guppy_args = {
        "expt_dir": expt_dir, 
        "config": config,
        "basecalling_method": basecalling_method,
        "barcoding_strategy": barcoding_strategy,
        "require_both_ends": "-b" if both_ends else "",
        "barcode_strict": "-s" if strict else "" 
    }

    # Formatting arguments for the rest of nomadic pipeline
    pipe_args = {"expt_dir": expt_dir, "config": config, "array_str": array_str}

    load_format_write(GUPPY_PIPELINE, output_dir / GUPPY_PIPELINE.name, **guppy_args)
    load_format_write(BAR_PIPELINE, output_dir / BAR_PIPELINE.name, **pipe_args)
    load_format_write(EXPT_PIPELINE, output_dir / EXPT_PIPELINE.name, **pipe_args)
    load_format_write(DWNSAMP_PIPELINE, output_dir / DWNSAMP_PIPELINE.name, **pipe_args)

    submission_path =  output_dir / SUBMIT_PIPELINE.name
    load_format_write(SUBMIT_PIPELINE, output_dir / SUBMIT_PIPELINE.name, run_dir=output_dir)
    submission_path.chmod(0o777)
    print(f"Final submission script written to: {output_dir / SUBMIT_PIPELINE.name}")

if __name__ == "__main__":
    main()

