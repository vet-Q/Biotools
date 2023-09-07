import click
from nomadic.lib.generic import produce_dir
from .basecalling import FLOW_CELLS, BASECALL_METHODS, run_guppy_basecaller
from .barcoding import BARCODING_KIT_MAPPING, run_guppy_barcode


# ================================================================
# Parameters
#
# ================================================================


ONLY_PASS = True  # only demultiplex .fastq that pass guppy quality control


# ================================================================
# Run basecalling with guppy
#
# ================================================================


@click.command(short_help="Basecall with guppy.")
@click.option(
    "-e",
    "--expt_dir",
    type=str,
    required=True,
    help="Path to experiment directory.",
)
@click.option(
    "-f",
    "--flow_cell",
    type=click.Choice(FLOW_CELLS),
    default="R10",
    help="Flow Cell chemistry.",
)
@click.option(
    "-m",
    "--basecalling_method",
    type=click.Choice(BASECALL_METHODS),
    default="hac",
    help="Basecalling method; super, high or fast accuracy.",
)
@click.option(
    "-q",
    "--min_qscore",
    type=int,
    default=12,
    help="Minimum Q-score for read to pass basecalling quality filter",
)
def basecall(expt_dir, flow_cell, basecalling_method, min_qscore):
    """
    Run guppy basecalling on an experiment

    """

    # Define directories
    input_dir = f"{expt_dir}/minknow"
    output_dir = produce_dir(expt_dir, "guppy", basecalling_method)

    # Print to stdout
    print(f"Experiment dir: {expt_dir}")
    print(f"Input dir: {input_dir}")
    print(f"Flow Cell: {flow_cell}")
    print(f"Basecalling method: {basecalling_method}")
    print(f"Minimum Q-score: {min_qscore}")
    print(f"Output dir: {output_dir}")

    # Run guppy basecalling
    run_guppy_basecaller(
        input_dir=input_dir, output_dir=output_dir, method=basecalling_method
    )


# ================================================================
# Run demultiplexing with guppy
#
# ================================================================


@click.command(short_help="Demultiplex with guppy.")
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
    type=click.Choice(BASECALL_METHODS),
    default="hac",
    help="Basecalling method; super, high or fast accuracy.",
)
@click.option(
    "-k",
    "--barcoding_strategy",
    type=click.Choice(BARCODING_KIT_MAPPING),
    default="R10_native96",
    help="Barcoding strategy, human-readable names that map to ONT kits.",
)
@click.option(
    "-b",
    "--both_ends",
    is_flag=True,
    default=False,
    help="Require both ends to have barcode?",
)
@click.option(
    "-s",
    "--strict",
    is_flag=True,
    default=False,
    help="Only classify barcodes if alignment score exceeds a strict threshold.",
)
def barcode(expt_dir, basecalling_method, barcoding_strategy, both_ends, strict):
    """
    Run guppy demultiplexing on .fastq files

    """

    # SELECT KIT
    barcode_kits = BARCODING_KIT_MAPPING[barcoding_strategy]

    # CREATE OUTPUT DIRECTORY
    input_dir = f"{expt_dir}/guppy/{basecalling_method}"
    if ONLY_PASS:
        fastq_input_dir = f"{input_dir}/pass"

    dir_name = "both_ends" if both_ends else "single_end"
    if strict:
        dir_name += "_strict"

    output_dir = produce_dir(input_dir, dir_name)

    # PRINT TO STDOUT
    print("Inputs")
    print(f"  Experiment dir.: {expt_dir}")
    print(f"  Basecalling method: {basecalling_method}")
    print(f"  Input directory: {input_dir}")
    if both_ends:
        print("  Requiring both ends to be barcoded.")
    else:
        print("  Requiring only one end to be barcoded.")
    if strict:
        print("  Implementing a strict barcode classification.")
    print(f"  Barcode kits: {barcode_kits}")
    print("  Will run with standard configuration.")
    print(f"  Output directory: {output_dir}")
    print("Done.")
    print("")

    # RUN GUPPY
    print("Running guppy barcoder...")
    run_guppy_barcode(
        fastq_input_dir=fastq_input_dir,
        barcode_kits=barcode_kits,
        output_dir=output_dir,
        both_ends=both_ends,
        strict=strict
    )
    print("Done.")
    print("")
