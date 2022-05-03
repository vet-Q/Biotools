
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

cli.add_command(search)

if __name__ == "__main__":
    cli()