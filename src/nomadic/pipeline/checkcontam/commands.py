import click
from nomadic.pipeline.cli import experiment_options, barcode_option

@click.command(short_help="Check contamination rate across an experiment.")
@experiment_options
def checkcontam(expt_dir, config):
    """
    Quickly call variants and annotate them with a given
    variant calling method

    """
    from .main import checkcontam

    checkcontam(expt_dir, config)

