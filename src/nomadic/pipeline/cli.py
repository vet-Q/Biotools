import click


# ================================================================
# Decorators for commonly used options
#
# ================================================================


def experiment_options(fn):
    """
    Wrapper for Click arguments used to specify the experiment,
    name -e <expt_dir> and -c <config_file>

    """
    fn = click.option(
        "-c",
        "--config",
        type=str,
        default="configs/default.ini",
        help="Path to NOMADIC configuration (.ini) file.",
    )(fn)
    fn = click.option(
        "-e",
        "--expt_dir",
        type=str,
        required=True,
        help="Path to experiment directory.",
    )(fn)
    return fn


def barcode_option(fn):
    """
    Wrapper for Click argument used to specify a specific
    barcode

    """
    fn = click.option(
        "-b",
        "--barcode",
        type=int,
        help="Optionally run command for only a single barcode.",
    )(fn)
    return fn


# ================================================================
# Entry point for all commands
#
# ================================================================


@click.group()
def cli():
    """
    NOMADIC: A pipeline for analysis of malaria long-read data

    """
    pass


# ================================================================
# Individual commands
#
# ================================================================


from .guppy.commands import barcode
from .map.commands import map
from .remap.commands import remap
from .qcbams.commands import qcbams
from .targets.commands import targets
from .calling.commands import call
from .find.commands import find

cli.add_command(barcode)
cli.add_command(map)
cli.add_command(remap)
cli.add_command(qcbams)
cli.add_command(targets)
cli.add_command(call)
cli.add_command(find)

from .bmrc.commands import bmrc

cli.add_command(bmrc)


if __name__ == "__main__":
    cli()
