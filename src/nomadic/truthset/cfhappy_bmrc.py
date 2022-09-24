# Idea is to create a script
# That creates submissions to BMRC
# Of truthset cfhappy
# Iterating over:

# Barcodes
# Methods

# Each barcode X method combination sent to a different
# node on the cluster

# Flexible enough to be extended to OGC comparison in near future
# How?
# - Mapping between barcodes and truth VCFs needs to be flexible
# - This is a simple dictionary
# - Store as a .json
# - With a small script to write the .json


import re
import click
import json
from itertools import product
from nomadic.pipeline.cli import experiment_options, barcode_option
from nomadic.pipeline.bmrc.commands import BmrcScriptBuilder


# ================================================================
# Main script
#
# ================================================================


@click.command(short_help="BMRC submission for hap.py.")
@experiment_options
@click.option(
    "-j",
    "--json_mapping",
    type=click.Path(),
    required=True,
    help="JSON file containing barcode to truth VCF mapping.",
)
@click.option(
    "-m",
    "--methods",
    type=str,
    required=True,
    help="Comma seperated list of calling methods.",
)
@click.option(
    "-f",
    "--bed_path",
    type=click.Path(),
    required=True,
    help="Path to BED file defining regions in which"
    " all variants in `truth_vcf` will be considered true positives.",
)
@click.option(
    "-s",
    "--stratification",
    type=click.Path(),
    required=False,
    default=None,
    help="Path to TSV file for stratifications.",
)
@click.option(
    "--downsample",
    is_flag=True,
    help="Compare against downsampled vcfs (e.g. in /downsample).",
)
def cfhappybmrc(
    expt_dir, config, json_mapping, methods, bed_path, stratification=None, downsample=False
):
    """
    Generate a submission script for hap.py comparisons

    """

    # Parse inputs
    barcode_mapping = json.load(open(json_mapping, "r"))
    method_list = [m.strip() for m in methods.split(",")]

    # Define submission file name
    submission_fn = "submit_happy-id.sh"

    # Write submission file
    with open(submission_fn, "w") as submit_file:
        for (barcode, truth_vcf), method in product(
            barcode_mapping.items(), method_list
        ):

            # Extract barcode integer
            barcode_int = int(re.match("^barcode([0-9]+)", barcode).group(1))

            # Construct command
            cmd = "truthset cfhappy"
            cmd += f" -e {expt_dir}"
            cmd += f" -c {config}"
            cmd += f" -m {method}"
            cmd += f" -b {barcode_int}"
            cmd += f" -t {truth_vcf}"
            cmd += f" -f {bed_path}"
            if stratification is not None:
                cmd += f" --stratification {stratification}"
            if downsample:
                cmd += " --downsample"

            # Write submission script
            builder = BmrcScriptBuilder(
                script=cmd,
                job_name=f"ch-{barcode_int:02d}-{method}"
            )
            builder.create_bmrc_script(kwargs=None)
            submit_file.write(builder.get_qsub_statement())