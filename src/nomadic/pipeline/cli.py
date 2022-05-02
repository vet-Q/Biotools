import click
from nomadic.pipeline.calling.callers import caller_collection
from nomadic.pipeline.submit_bmrc import PIPELINE_PATH
from nomadic.pipeline.guppy_barcode import BARCODING_KIT_MAPPING


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
    

# def common_options(fn):
#     fn = click.option(
#         "-b",
#         "--barcode",
#         type=int,
#         help="Optionally run command for only a single barcode.",
#     )(fn)
#     fn = click.option(
#         "-c",
#         "--config",
#         type=str,
#         default="configs/default.ini",
#         help="Path to NOMADIC configuration (.ini) file.",
#     )(fn)
#     fn = click.option(
#         "-e",
#         "--expt_dir",
#         type=str,
#         required=True,
#         help="Path to experiment directory.",
#     )(fn)
#     return fn


# ================================================================
# Entry point for all commands
#
# ================================================================


# I think this is where I set logging verbosity
# Then we want that verbosity level to pass to the sub-modules, probably
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


@cli.command(short_help="Run complete pipeline.")
@experiment_options
@barcode_option
def runall(expt_dir, config, barcode):
    """
    Run the complete NOMADIC pipeline

    """
    from nomadic.pipeline import (
        map_pf,
        remap_to_hs,
        qc_bams,
        qc_bams_overview,
        target_extraction,
        target_extraction_overview,
    )

    map_pf.main(expt_dir, config, barcode)
    remap_to_hs.main(expt_dir, config, barcode)
    qc_bams.main(expt_dir, config, barcode)
    qc_bams_overview.main(expt_dir, config)
    target_extraction.main(expt_dir, config, barcode)
    target_extraction_overview.main(expt_dir, config)


@cli.command(short_help="Demultiplex with guppy.")
@click.option(
        "-e",
        "--expt_dir",
        type=str,
        required=True,
        help="Path to experiment directory.",
    )
@click.option(
    "-m",
    "--basecalling_method",
    type=click.Choice(["hac", "fast"]),
    default="hac",
    help="Basecalling method, high accuracy or fast."
)
@click.option(
    "-k",
    "--barcoding_strategy",
    type=click.Choice(BARCODING_KIT_MAPPING),
    default="native24",
    help="Barcoding strategy, human-readable names that map to ONT kits."
)
@click.option(
    "-b", 
    "--both_ends",
    is_flag=True,
    default=False,
    help="Require both ends to have barcode?"
)
def barcode(expt_dir, basecalling_method, barcoding_strategy, both_ends):
    from nomadic.pipeline import guppy_barcode
    guppy_barcode.main(expt_dir, basecalling_method, barcoding_strategy, both_ends)


@cli.command(short_help="Map to P.f. reference.")
@experiment_options
@barcode_option
def map(expt_dir, config, barcode):
    """
    Map .fastq files to the P. falciparum reference genome

    """
    from nomadic.pipeline import map_pf

    map_pf.main(expt_dir, config, barcode)


@cli.command(short_help="Map unmapped reads to H.s.")
@experiment_options
@barcode_option
def remap(expt_dir, config, barcode):
    """
    Some reads may fail to map to the P.f. reference genome after
    running `nomadic map`

    `nomadic remap` tries to map these reads to the H.s. reference
    genome

    """
    from nomadic.pipeline import remap_to_hs

    remap_to_hs.main(expt_dir, config, barcode)


@cli.command(short_help="QC analysis of .bam files.")
@experiment_options
@barcode_option
@click.option(
    "--overview", is_flag=True, help="Produce an overview across all barcodes."
)
def qcbams(expt_dir, config, barcode, overview):
    """
    Run a quality control analysis of .bam files generated from
    `nomadic map` and `nomadic remap`

    """
    if overview:
        from nomadic.pipeline import qc_bams_overview

        qc_bams_overview.main(expt_dir, config)
    else:
        from nomadic.pipeline import qc_bams

        qc_bams.main(expt_dir, config, barcode)


@cli.command(short_help="Analyse amplicon targets.")
@experiment_options
@barcode_option
@click.option(
    "--overview", is_flag=True, help="Produce an overview across all barcodes."
)
def targets(expt_dir, config, barcode, overview):
    """
    Run analyses of a specific set of amplicon targets

    TODO:
    - Probably want an additional argument here specifying targets
    - Would add flexibility

    """
    if overview:
        from nomadic.pipeline import target_extraction_overview

        target_extraction_overview.main(expt_dir, config)
    else:
        from nomadic.pipeline import target_extraction

        target_extraction.main(expt_dir, config, barcode)


@cli.command(short_help="Call variants across targets.")
@experiment_options
@barcode_option
@click.option(
    "-m",
    "--method",
    type=click.Choice(caller_collection),
    required=True,
    help="Variant calling method to use.",
)
@click.option(
    "--downsample", is_flag=True, help="Produce an overview across all barcodes."
)
def call(expt_dir, config, barcode, method, downsample):
    """
    Run analyses of a specific set of amplicon targets

    TODO:
    - Probably want an additional argument here specifying targets
    - Would add flexibility

    """
    from nomadic.pipeline.calling import commands
    if downsample:
        pass
    else:
        commands.call(expt_dir, config, barcode, method)


@cli.command(short_help="Find mutations of interest.")
@experiment_options
@barcode_option
@click.option(
    "-m",
    "--method",
    type=click.Choice(caller_collection),
    required=True,
    help="Variant calling method to use.",
)
def find(expt_dir, config, barcode, method):
    """
    Find a set of mutations in your amplicon set


    """
    from nomadic.pipeline.find import commands
    commands.main(expt_dir, config, barcode, method)


@cli.command(short_help="Build BMRC pipeline submission.")
@experiment_options
@click.option(
    "-p", "--pipeline", 
    type=str,
    default=PIPELINE_PATH,
    help="Path to BMRC pipeline (.ini) file.")
def bmrc(expt_dir, config, pipeline):
    """
    Build necessary submission scripts to run pipeline
    on the BMRC cluster

    """
    from nomadic.pipeline import submit_bmrc

    submit_bmrc.main(expt_dir, config, pipeline)


from nomadic.pipeline.find.commands import find

cli.add_command(find)


if __name__ == "__main__":
    cli()
