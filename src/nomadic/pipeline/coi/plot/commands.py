import click
from ..trim.targets import TARGET_COLLECTION
from .main import main
from nomadic.pipeline.cli import experiment_options, barcode_option


@click.command(short_help="Plot COI results.")
@experiment_options
@barcode_option
@click.option(
    "-t",
    "--target_gene",
    type=click.Choice(TARGET_COLLECTION),
    default="MSP2",
    help="Trim mapped reads to this target gene."
)
def plot(expt_dir, config, barcode, target_gene):
    """
    Plot results of COI analyses.

    Note: must run `trim`, `panmap`, `overlap`, and 
    `align` first.
    
    """
    main(expt_dir, config, barcode, target_gene)