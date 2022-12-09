import click
from ..trim.targets import TARGET_COLLECTION
from .main import main
from nomadic.pipeline.cli import experiment_options, barcode_option


@click.command(short_help="Map to a panel of Pf strains.")
@experiment_options
@barcode_option
@click.option(
    "-t",
    "--target_gene",
    type=click.Choice(TARGET_COLLECTION),
    default="MSP2",
    help="Target gene reads to map."
)
def panmap(expt_dir, config, barcode, target_gene):
    """
    Map reads from a `target_gene` to a panel of P.f.
    referennce strains
    
    """
    main(expt_dir, config, barcode, target_gene)