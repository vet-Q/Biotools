
import click

# ================================================================
# Entry point
#
# ================================================================


@click.group()
def cli():
    """
    TRUTHSET: Define true variants across malaria strains

    """
    pass


# ================================================================
# Individual commands
#
# ================================================================


from .search import search
from .fasta import fasta
from .msacall import msacall

cli.add_command(search)
cli.add_command(fasta)
cli.add_command(msacall)

if __name__ == "__main__":
    cli()