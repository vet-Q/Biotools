import os
import re
import subprocess
import allel
import pandas as pd
import numpy as np

from dataclasses import dataclass
from typing import Dict, List
from nomadic.lib.generic import produce_dir
from nomadic.lib.parsing import build_parameter_dict


# --------------------------------------------------------------------------------
# Load VCFs
#
# --------------------------------------------------------------------------------


def load_cfhappy_vcf_as_dict(vcf_path: str) -> Dict:
    """
    Load cfhappy output VCF file as a dictionary

    """
    vcf = allel.read_vcf(
        vcf_path,
        fields=[
            "samples",
            "variants/CHROM",
            "variants/POS",
            "variants/REF",
            "variants/ALT",
            "variants/BS",  # Benchmark superlocus
            "calldata/BD",
        ],
    )
    return vcf


def load_cfhappy_vcf_as_dataframe(vcf_path: str) -> pd.DataFrame:
    """
    Load cfhappy output VCF file as a dataframe

    """
    # Load as dictionary
    vcf = load_cfhappy_vcf_as_dict(vcf_path)

    # Convert to dataframe
    vcf_df = pd.DataFrame(
        {
            "CHROM": vcf["variants/CHROM"],
            "POS": vcf["variants/POS"],
            "REF": vcf["variants/REF"],
            "ALT": vcf["variants/ALT"][:, 0],
            "BS": vcf["variants/BS"],
            "TRUTH_CALL": vcf["calldata/BD"][:, 0],
            "QUERY_CALL": vcf["calldata/BD"][:, 1],
        }
    )

    return vcf_df


def convert_cfhappy_call(truth_call: str, query_call: str):
    """
    Convert a VCF call from cfhappy to binary

    """

    # Store like so
    @dataclass
    class CallingResult:
        TP: float
        FP: float
        FN: float

    # Classify
    if truth_call == "TP" and query_call == "TP":
        return CallingResult(1, 0, 0)
    elif truth_call == "." and query_call == "FP":
        return CallingResult(0, 1, 0)
    elif truth_call == "FN" and query_call == ".":
        return CallingResult(0, 0, 1)
    elif truth_call == "FN" and query_call == "FP":  # het instead of hom
        return CallingResult(0.5, 0, 0.5)
    else:
        raise ValueError(
            f"Unable to classify, TRUTH: {truth_call}, QUERY: {query_call}"
        )


def add_vcf_call_columns(vcf_df):
    """
    Add binary cfhappy calling columns

    """

    calls_df = pd.DataFrame(
        [
            convert_cfhappy_call(t, q)
            for t, q in zip(vcf_df["TRUTH_CALL"], vcf_df["QUERY_CALL"])
        ]
    )

    assert (calls_df.sum(1) == 1).all()

    return pd.concat([vcf_df, calls_df], axis=1)


# --------------------------------------------------------------------------------
# Main script
#
# --------------------------------------------------------------------------------


def main(expt: str, config: str, barcodes: List[str], method: str) -> None:
    """
    Aggregate cfhappy VCF information for a given `expt`, `config`,
    set of `barcodes` and variant calling `method`

    """

    # Load experiment parameters
    params = build_parameter_dict(expt, config, barcode=None)

    vcf_dfs = []
    for barcode in barcodes:
        print(f"Aggregating VCF information for {barcode}...")

        # Create directories
        cfhappy_dir = f"{params['barcodes_dir']}/{barcode}/cfhappy/{method}"
        output_dir = produce_dir(cfhappy_dir, "positional_analysis")

        # Load cfhappy VCFs
        all_targets_pattern = (
            f"{barcode}.n([0-9]{'{4}'}).r([0-9]{'{3}'}).all_targets.sorted.gz.vcf.gz"
        )

        # Iterate over files in `cfhappy_dir`...
        for f in os.listdir(cfhappy_dir):
            # COLLECT IF VCF
            match_result = re.match(all_targets_pattern, f)
            if match_result is None:
                continue

            # REDUCE TO SNPs in CONFIDENT REGIONS
            # File names
            vcf_path = f"{cfhappy_dir}/{f}"
            replace_str = ".sorted.gz.vcf.gz"
            assert vcf_path.endswith(replace_str)
            clean_vcf_path = f"{output_dir}/{f.replace(replace_str, '.cleaned.vcf.gz')}"

            # Run command
            cmd = f"bcftools view {vcf_path}"
            cmd += " --type='snps'"
            cmd += " -e 'BD=\"UNK\"'"
            cmd += " -Ob"
            cmd += f" -o {clean_vcf_path}"
            # print(cmd)
            subprocess.run(cmd, check=True, shell=True)

            # LOAD VCF
            vcf_df = load_cfhappy_vcf_as_dataframe(clean_vcf_path)
            # Compute TP, FP, FN rates
            vcf_df = add_vcf_call_columns(vcf_df)

            # Annotate
            vcf_df.insert(0, "barcode", barcode)
            vcf_df.insert(1, "method", method)
            vcf_df.insert(2, "n_reads", int(match_result.group(1)))
            vcf_df.insert(3, "rep", int(match_result.group(2)))

            # STORE
            vcf_dfs.append(vcf_df)

    # Create final output directory
    assert len(vcf_dfs) == (80 * len(barcodes))
    combined_vcf_df = pd.concat(vcf_dfs)

    # Write
    combined_vcf_df.to_csv(
        f"summary_tables/fig4_aggregated-vcfs-{method}.csv", index=False
    )


if __name__ == "__main__":
    expt = "experiments/2021-11-14_strain-validation-flongle-lfb"
    config = "configs/guppy-hac-se.ini"
    barcodes = ["barcode02", "barcode03", "barcode04"]
    method = "clair3sing"
    main(expt, config, barcodes, method)
