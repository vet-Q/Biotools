import click
from nomadic.pipeline.cli import experiment_options, barcode_option

# Defaults
AMPLICON_BED_PATH = (
    "resources/truthsets/stratifications/multiplex.03.greedy.confident_amplicons.bed"
)
CDS_BED_PATH = "resources/truthsets/stratifications/multiplex.03.greedy.cds.bed "

# Command
@click.command(short_help="Analyse read depth profile.")
@experiment_options
@barcode_option
@click.option(
    "-a",
    "--amplicon_bed_path",
    type=str,
    default=AMPLICON_BED_PATH,
    help="Path to BED file defining amplicon regions.",
)
@click.option(
    "-d",
    "--cds_bed_path",
    type=str,
    default=CDS_BED_PATH,
    help="Path to BED file defining amplicon coding sequences regions.",
)
def depth(expt_dir, config, barcode, amplicon_bed_path, cds_bed_path):
    """
    Analyse depth profiles across a set of amplicons

    """
    from .main import depth

    depth(expt_dir, config, barcode, amplicon_bed_path, cds_bed_path)
