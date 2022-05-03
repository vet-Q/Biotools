
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

cli.add_command(search)
cli.add_command(fasta)

if __name__ == "__main__":
    cli()