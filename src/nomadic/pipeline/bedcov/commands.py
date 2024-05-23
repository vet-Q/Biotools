import click
from nomadic.pipeline.cli import experiment_options, barcode_option


@click.command(short_help="Compute coverage across regions.")
@experiment_options
@barcode_option
@click.option(
    "-r", "--bed_path", type=str, help="Path to BED file for computing coverage."
)
@click.option(
    "--overview", is_flag=True, help="Produce an overview across all barcodes."
)
def bedcov(expt_dir, config, barcode, bed_path, overview):
    """
    Compute coverage across a set of regions defined by a BED file

    """

    from .main import bedcov

    bedcov(expt_dir, config, bed_path, barcode, overview)
