import os
import click

from nomadic.lib.generic import produce_dir
from nomadic.lib.parsing import build_parameter_dict
from nomadic.pipeline.cli import experiment_options, barcode_option


PIPELINE = "slurm/templates/sbatch-pipeline.sh"
STEP_DIR = "slurm/templates/steps"
RUNS_DIR = "slurm/runs"


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

    # Get individual steps
    step_scripts = [
        os.path.join(f"{STEP_DIR}/{s}") 
        for s in os.listdir(STEP_DIR)
        if s.endswith(".sh") and s.startswith("nomadic")
    ]
    print(f"Found {len(step_scripts)} to add to pipeline...")

    # Write individual steps
    for script in step_scripts:
        # Format
        script_str = "".join(open(script, "r").readlines())
        script_str_formatted = script_str.format(
            expt_dir=expt_dir,
            config=config,
            array_str=array_str
        )

        # Write
        script_name = os.path.basename(script)
        with open(f"{output_dir}/{script_name}", "w") as output_script:
            output_script.write(script_str_formatted)

    # Write pipeline
    pipeline_str = "".join(open(PIPELINE, "r").readlines())
    pipeline_str_formatted = pipeline_str.format(
        run_dir=output_dir
    )
    with open(f"{output_dir}/{os.path.basename(PIPELINE)}", "w") as output_pipeline:
        output_pipeline.write(pipeline_str_formatted)


if __name__ == "__main__":
    main()