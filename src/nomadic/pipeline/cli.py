import click


# ================================================================
# Decorators for commonly used options
#
# ================================================================


def common_options(fn):
    fn = click.option(
        "-b",
        "--barcode",
        type=int,
        help="Optionally run command for only a single barcode.",
    )(fn)
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
@common_options
def runall(expt_dir, config, barcode):
    """
    Run the complete NOMADIC pipeline

    """
    from nomadic.pipeline import map_pf, remap_to_hs, qc_bams, target_extraction

    map_pf.main(expt_dir, config, barcode)
    remap_to_hs.main(expt_dir, config, barcode)
    qc_bams.main(expt_dir, config, barcode)
    target_extraction.main(expt_dir, config, barcode)


@cli.command(short_help="Map to P.f. reference.")
@common_options
def map(expt_dir, config, barcode):
    """
    Map .fastq files to the P. falciparum reference genome

    """
    from nomadic.pipeline import map_pf

    map_pf.main(expt_dir, config, barcode)


@cli.command(short_help="Map unmapped reads to H.s.")
@common_options
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
@common_options
@click.option("--overview", is_flag=True, help="Produce an overview across all barcodes.")
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
@common_options
@click.option("--overview", is_flag=True, help="Produce an overview across all barcodes.")
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


@cli.command(short_help="Build BMRC pipeline submission.")
@common_options
def bmrc(expt_dir, config, barcode):
    """
    Build necessary submission scripts to run pipeline
    on the BMRC cluster

    """
    from nomadic.pipeline import submit_bmrc

    submit_bmrc.main(expt_dir, config, barcode)


if __name__ == "__main__":
    cli()
