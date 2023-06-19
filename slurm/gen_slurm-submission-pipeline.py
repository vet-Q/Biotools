import os
import click

from nomadic.lib.generic import produce_dir
from nomadic.lib.parsing import build_parameter_dict
from nomadic.pipeline.cli import experiment_options, barcode_option


BAR_PIPELINE = "slurm/templates/pipelines/nomadic-barcode.sh"
EXPT_PIPELINE = "slurm/templates/pipelines/nomadic-experiment.sh"
SUBMIT_PIPELINE = "slurm/templates/nomadic-submit.sh"
RUNS_DIR = "slurm/runs"


def load_format_write(input_file, output_file, **kwargs):
    """Load a file that has named formating fields, e.g. {job_name}, format it, and write"""
    
    input_str = "".join(open(input_file, "r").readlines())
    output_str = input_str.format(**kwargs)
    with open(output_file, "w") as output:
        output.write(output_str)



@click.command(short_help="Generate a NOMADIC submission script for Raven.")
@experiment_options
@barcode_option
def main(expt_dir: str, config: str, barcode: str):
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
    expt_name = os.path.basename(expt_dir)
    output_dir = produce_dir(RUNS_DIR, expt_name)

    # Format
    pipe_args = {"expt_dir": expt_dir, "config": config, "array_str": array_str}
    load_format_write(BAR_PIPELINE, f"{output_dir}/{os.path.basename(BAR_PIPELINE)}", **pipe_args)
    load_format_write(EXPT_PIPELINE, f"{output_dir}/{os.path.basename(EXPT_PIPELINE)}", **pipe_args)
    load_format_write(SUBMIT_PIPELINE, f"{output_dir}/{os.path.basename(SUBMIT_PIPELINE)}", run_dir=output_dir)


if __name__ == "__main__":
    main()

