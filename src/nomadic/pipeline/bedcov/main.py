import os
import warnings
import pandas as pd
from nomadic.lib.process_bams import samtools_bedcov
from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.references import PlasmodiumFalciparum3D7


# ================================================================================ #
# Individual Barcode
#
# ================================================================================ #


def bedcov_single(
    expt_dir: str, config: str, bed_path: str, barcode: str = None
) -> None:
    """
    Compute coverage across a BED file for individual barcodes

    """
    # PARSE INPUTS
    script_descrip = "NOMADIC: Compute coverage across regions"
    t0 = print_header(script_descrip)
    script_dir = "bedcov"  # this will be easier to iterate overs
    params = build_parameter_dict(expt_dir, config, barcode)

    # Define reference genome
    reference = PlasmodiumFalciparum3D7()

    # Focus on a single barcode, if specified
    if "focus_barcode" in params:
        params["barcodes"] = [params["focus_barcode"]]

    # Iterate over barcodes
    for barcode in params["barcodes"]:
        # Define input and output directory
        print("-" * 80)
        print(f"Running for: {barcode}")
        barcode_dir = f"{params['barcodes_dir']}/{barcode}"
        input_dir = f"{barcode_dir}/bams"
        output_dir = produce_dir(barcode_dir, script_dir)

        # Path to *complete* bam file
        bam_path = f"{input_dir}/{barcode}.{reference.name}.final.sorted.bam"
        csv_path = f"{output_dir}/{barcode}.{reference.name}.bedcov.csv"

        # Compute BED coverage
        # -> Does not add any sample or barcode information natively
        samtools_bedcov(bam_path=bam_path, bed_path=bed_path, output_csv=csv_path)

    print_footer(t0)


# ================================================================================ #
# Merge across barcodes
#
# ================================================================================ #


def bedcov_merge(expt_dir: str, config: str, bed_path: str) -> None:
    """
    Merge results from `bedcov_single`

    - Annotate the CSV files with barcode name
    - Could also compute two summary tables
        -> Each barcode, across amplicons
        -> Each amplicon, across barcodes

    """
    # PARSE INPUTS
    script_descrip = "NOMADIC: Compute coverage across regions"
    t0 = print_header(script_descrip)
    params = build_parameter_dict(expt_dir, config)
    script_dir = "bedcov"  # this will be easier to iterate overs
    output_dir = produce_dir(params["nomadic_dir"], script_dir)

    # Define reference genome
    reference = PlasmodiumFalciparum3D7()

    # Iterate over barcodes
    dfs = []
    for barcode in params["barcodes"]:
        # Define input and output directory
        print("-" * 80)
        print(f"Running for: {barcode}")
        barcode_dir = f"{params['barcodes_dir']}/{barcode}"
        input_dir = f"{barcode_dir}/bedcov"

        # Path to DF
        csv_path = f"{input_dir}/{barcode}.{reference.name}.bedcov.csv"

        # Skip if missing
        if not os.path.exists(csv_path):
            warnings.warn(f"No CSV file found at {csv_path}! Skipping.")
            continue

        # Load, annotate, append
        barcode_df = pd.read_csv(csv_path)
        barcode_df.insert(0, "barcode", barcode)
        dfs.append(barcode_df)

    # Concatenate
    print("Merging all coverage data...")
    final_df = pd.concat(dfs)
    print(f"  Found {len(final_df['name'].unique())} amplicons.")
    print(f"  Found {len(final_df['barcode'].unique())} barcodes.")
    print(f"  Resulting in information for {final_df.shape[0]} sets of (amplicon, sample).")

    # Write
    output_csv = f"{output_dir}/summary.bedcov.csv"
    print(f"  Writing to: {output_csv}")
    final_df.to_csv(output_csv)
    print("Done.")

    print_footer(t0)


# ================================================================================ #
# Interface
#
# ================================================================================ #


def bedcov(
    expt_dir: str, config: str, barcode: str, bed_path: str, overview: bool = False
) -> None:
    if overview:
        bedcov_merge(expt_dir, config, bed_path)
    else:
        bedcov_single(expt_dir, config, barcode, bed_path)
