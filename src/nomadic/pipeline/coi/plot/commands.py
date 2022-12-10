import click
from ..trim.targets import TARGET_COLLECTION
from .main import main
from .overview import plot_overview
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
@click.option(
    "--overview", is_flag=True, help="Produce an overview across all barcodes."
)
def plot(expt_dir, config, barcode, target_gene, overview):
    """
    Plot results of COI analyses.

    Note: must run `trim`, `panmap`, `overlap`, and 
    `align` first.
    
    """

    if overview:
        plot_overview(expt_dir, config, target_gene)
    else:
        main(expt_dir, config, barcode, target_gene)