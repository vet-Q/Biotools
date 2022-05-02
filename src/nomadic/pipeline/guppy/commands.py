import click
from nomadic.lib.generic import produce_dir
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


# ADD CLICK
def basecall():
    pass


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
    type=click.Choice(["hac", "fast"]),
    default="hac",
    help="Basecalling method, high accuracy or fast.",
)
@click.option(
    "-k",
    "--barcoding_strategy",
    type=click.Choice(BARCODING_KIT_MAPPING),
    default="native24",
    help="Barcoding strategy, human-readable names that map to ONT kits.",
)
@click.option(
    "-b",
    "--both_ends",
    is_flag=True,
    default=False,
    help="Require both ends to have barcode?",
)
def barcode(expt_dir, basecalling_method, barcoding_strategy, both_ends):
    """
    Run guppy demultiplexing on .fastq files

    """

    # SELECT KIT
    barcode_kits = BARCODING_KIT_MAPPING[barcoding_strategy]

    # CREATE OUTPUT DIRECTORY
    input_dir = f"{expt_dir}/guppy/{basecalling_method}"
    if ONLY_PASS:
        fastq_input_dir = f"{input_dir}/pass"
    output_dir = produce_dir(input_dir, "both_ends" if both_ends else "single_end")

    # PRINT TO STDOUT
    print("Inputs")
    print(f"  Experiment dir.: {expt_dir}")
    print(f"  Basecalling method: {basecalling_method}")
    print(f"  Input directory: {input_dir}")
    if both_ends:
        print(f"  Requiring both ends to be barcoded.")
    else:
        print("  Requiring only one end to be barcoded.")
    print(f"  Barcode kits: {barcode_kits}")
    print(f"  Will run with standard configuration.")
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
    )
    print("Done.")
    print("")
