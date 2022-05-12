
import click


# TODO:
# - Advatages, easy to do add to the truthset from minimap2, or from other illumina data
# - Can probably run quite a similar pipeline for the Illumina comparison data:
#   



# Additional command:
#   truthset merge -e <expt_dir> -c <config> -m <method> --truthset <path>.vcf --downsample
#
#   Merge all of the VCFs for all barcodes for a particular experiment
#
#   if --downsample, then merge material in downsample directory
#   Need to come up with good standardised naming for the samples
#   We will also need to restrict to only SNPs, i.e. we want to remove structural variants
#
#   if --truthset <path.vcf>, also merge in a specificed truthset <.vcf>
#
#
#   - Might want to restrict to specific barcodes, e.g. -b 1 -b 2 -b 3 -b 4, &c.
#   Output:  /calling/method/<downsample>/merged.<info>.vcf

# Update:
#  truthset msacall
#
#  We want this to output as a joint VCF now

# Final data-munging command:
#   truthset tabulate -p <path_to_merged.vcf>
#
#   This creates a tall data frame from the merged VCF. This tall data frame can be used to create all of the
#   necessary plots.

# Finally, I will want a "truthset run" command, that runs the full pipeline;
# Will make it easier to remember how things were generated later on.
#
#






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
from .merge import merge

cli.add_command(search)
cli.add_command(fasta)
cli.add_command(msacall)
cli.add_command(merge)

if __name__ == "__main__":
    cli()