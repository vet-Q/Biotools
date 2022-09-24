# Create BMRC submissions for cfhappy
# 2022/09/24, J.Hendry
#
# TODO:
# - Don't really want bmrc/bmrc_* for this; just the submission
# - So, would want long command line argument


import re
import os
import click
import json
import uuid
from itertools import product
from nomadic.pipeline.cli import experiment_options
from nomadic.truthset.cfhappy import happy_callers


# ================================================================
# Parameters
#
# ================================================================


HAPPY_BMRC = "bmrc/configs/bmrc_cfhappy-template.sh"
SUBMISSION_ID = str(uuid.uuid4())[:5]


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
    "-h",
    "--happy_caller",
    type=click.Choice(happy_callers),
    required=False,
    default="docker",
    help="Run hap.py by Docker or Singularity.",
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
    expt_dir,
    config,
    json_mapping,
    methods,
    bed_path,
    happy_caller,
    stratification=None,
    downsample=False,
):
    """
    Generate a submission script for hap.py comparisons

    """
    # Parse inputs
    print("Parsing inputs...")
    barcode_mapping = json.load(open(json_mapping, "r"))
    method_list = [m.strip() for m in methods.split(",")]

    # Define submission file name
    submission_fn = f"submit_happy-{SUBMISSION_ID}.sh"

    # Write submission file
    print("Writting submission file...")
    with open(submission_fn, "w") as submit_file:
        for (barcode, truth_vcf), method in product(
            barcode_mapping.items(), method_list
        ):

            # Extract barcode integer
            barcode_int = int(re.match("^barcode([0-9]+)", barcode).group(1))

            # Construct command
            cmd = f"qsub {HAPPY_BMRC}"
            cmd += f" -e {expt_dir}"
            cmd += f" -c {config}"
            cmd += f" -m {method}"
            cmd += f" -b {barcode_int}"
            cmd += f" -t {truth_vcf}"
            cmd += f" -f {bed_path}"
            cmd += f" -h {happy_caller}"
            if stratification is not None:
                cmd += f" --stratification {stratification}"
            if downsample:
                cmd += " --downsample"
            submit_file.write(f"{cmd}\n")
    os.chmod(submission_fn, 0o777)
    print(f" to: {submission_fn}")
    print("Done.")
