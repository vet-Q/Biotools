import click
from ..trim.targets import TARGET_COLLECTION
from .main import main
from .aligners import ALIGNER_COLLECTION
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
@click.option(
    "-m",
    "--max_reads",
    type=int,
    default=600,
    show_default=True,
    help="Maximum number of reads to pairwise align."
)
@click.option(
    "-a",
    "--algorithm",
    type=click.Choice(ALIGNER_COLLECTION),
    default="needleman_numba_banded_qscores",
    show_default=True,
    help="Pairwise alignment algorithm."
)
def align(expt_dir, config, barcode, target_gene, max_reads, algorithm):
    """
    Perform pairwise alignments for a collection
    of trimmed reads
    
    """
    main(expt_dir, config, barcode, target_gene, max_reads, algorithm)




    

    