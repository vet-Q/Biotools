import click
from ..trim.targets import TARGET_COLLECTION
from .main import main
from nomadic.pipeline.cli import experiment_options, barcode_option


@click.command(short_help="Look for overlaps between reads.")
@experiment_options
@barcode_option
@click.option(
    "-t",
    "--target_gene",
    type=click.Choice(TARGET_COLLECTION),
    default="MSP2",
    help="Target gene reads search."
)
def overlap(expt_dir, config, barcode, target_gene):
    """
    Use `minimap2` to look for overlaps between a set
    of reads deriving from a single `fastq`
    
    """
    main(expt_dir, config, barcode, target_gene)