import click
from .targets import TARGET_COLLECTION
from .main import main
from nomadic.pipeline.cli import experiment_options, barcode_option


@click.command(short_help="Trim to high-quality reads overlappng a target.")
@experiment_options
@barcode_option
@click.option(
    "-t",
    "--target_gene",
    type=click.Choice(TARGET_COLLECTION),
    default="MSP2",
    help="Trim mapped reads to this target gene."
)
def trim(expt_dir, config, barcode, target_gene):
    """
    Filter and trim all reads in a BAM file to overlap a `target_gene`, 
    then convert to FASTQ
    
    """
    main(expt_dir, config, barcode, target_gene)