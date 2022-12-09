import click
from ..trim.targets import TARGET_COLLECTION
from .main import main
from nomadic.pipeline.cli import experiment_options, barcode_option


@click.command(short_help="Perform pairwise alignment between reads.")
@experiment_options
@barcode_option
@click.option(
    "-t",
    "--target_gene",
    type=click.Choice(TARGET_COLLECTION),
    default="MSP2",
    help="Target gene reads to map."
)
def align(expt_dir, config, barcode, target_gene):
    """
    Perform pairwise alignments for a collection
    of trimmed reads
    
    """
    main(expt_dir, config, barcode, target_gene)




    

    