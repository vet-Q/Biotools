import click
from nomadic.pipeline.cli import cli
from .trim.commands import trim


@cli.group()
def coi():
    """
    COI estimation tools.

    """
    pass

coi.add_command(trim)