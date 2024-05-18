import click
from nomadic.pipeline.cli import experiment_options, barcode_option
from .callers import caller_collection


@click.command(short_help="Call variants across targets.")
@experiment_options
@barcode_option
@click.option(
        "-r",
        "--bed_path",
        type=str,
        help="Path to BED file for annotating amplicons."
    )
@click.option(
    "-m",
    "--method",
    type=click.Choice(caller_collection),
    default="bcftools",
    show_default=True,
    required=False,
    help="Variant calling method.",
)
def quickcall(expt_dir, config, barcode, bed_path, method):
    """
    Quickly call variants and annotate them with a given
    variant calling method

    """
    from .main import quickcall

    quickcall(expt_dir, config, barcode, bed_path, method)

