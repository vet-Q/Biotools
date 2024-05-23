import os
import warnings
import pandas as pd
import numpy as np

from typing import Dict
from functools import partial

from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.parsing import build_parameter_dict


# SETTINGS
NTC_INDICATOR = "NTC"  # anywhere within
NTC_COLUMN = "is_negative"


def create_ntc_column(metadata: pd.DataFrame, ntc_indicator: str) -> pd.DataFrame:
    """
    Create a boolean column indicating the presence of negative samples

    """
    return [ntc_indicator in sample_id for sample_id in metadata["sample_id"]]


def checkcontam(expt_dir: str, config: str):
    """
    Check for rates of contamination across an experiment
    where negative controls have been included

    """

    # PARSE INPUTS
    script_descrip = "NOMADIC: Assess the rate of contamination across an experiment"
    t0 = print_header(script_descrip)
    params = build_parameter_dict(expt_dir, config)
    script_dir = "contam"  # this will be easier to iterate overs
    input_dir = params["nomadic_dir"] + "/bedcov"
    output_dir = produce_dir(params["nomadic_dir"], script_dir)

    # Load the bed coverage dataframe
    print("Loading bed coverage data...")
    csv_path = f"{input_dir}/summary.bedcov.csv"
    bedcov_df = pd.read_csv(csv_path)

    # Load
    print("Loading metadata...")
    metadata = params["metadata"]

    # Define negative controls
    if NTC_COLUMN in metadata.columns:
        print(f"Found negative control column `{NTC_COLUMN}`.")
    else:
        print(f"Did not find negative control column `{NTC_COLUMN}`.")
        print(f"Inferring negative controls containing substring `{NTC_INDICATOR}`...")
        metadata.insert(
            metadata.shape[1], NTC_COLUMN, create_ntc_column(metadata, NTC_INDICATOR)
        )

    n_neg = metadata[NTC_COLUMN].sum()
    print(f"Found {n_neg} negative controls.")
    if n_neg == 0:
        raise ValueError("No negative controls -- cannot examine contamination levels.")

    # Merge negative control information
    n = bedcov_df.shape[0]
    merged_df = pd.merge(
        metadata[["barcode", "sample_id", "is_negative"]], bedcov_df, on="barcode"
    )
    assert n == merged_df.shape[0], "Lost samples during merging."

    # Compute NTC summary statistics
    total_cov = merged_df.query("not is_negative").groupby("name")["mean_cov"].sum()
    total_cov.name = "total_nonntc_cov"

    ntc_df = merged_df.query("is_negative")

    ntc_summary_df = (
        ntc_df.groupby("name")
        .agg(
            n_ntcs=pd.NamedAgg("mean_cov", len),
            mean_ntc_cov=pd.NamedAgg("mean_cov", np.mean),
            max_ntc_cov=pd.NamedAgg("mean_cov", np.max),
            min_ntc_cov=pd.NamedAgg("mean_cov", np.min),
        )
        .reset_index()
    )

    ntc_summary_df = pd.merge(left=ntc_summary_df, right=total_cov, on="name")

    ntc_summary_df.insert(
        ntc_summary_df.shape[1],
        "per_ntc_cov",
        100 * ntc_summary_df["mean_ntc_cov"] / ntc_summary_df["total_nonntc_cov"],
    )

    # Classify samples as PASS / FAIL
    ntc_summary_df.index = ntc_summary_df.name

    mean_ntc_dt = dict(ntc_summary_df["mean_ntc_cov"])

    def calc_per_mean_cov_ntc(row: pd.Series, mean_ntc_dt: Dict) -> float:
        return 100 * mean_ntc_dt[row["name"]] / (row["mean_cov"] + 0.01)

    def fold_ntc(row: pd.Series, mean_ntc_dt: Dict) -> float:
        return row["mean_cov"] / mean_ntc_dt[row["name"]]

    summary_funcs = {
        "per_mean_cov_ntc": partial(calc_per_mean_cov_ntc, mean_ntc_dt=mean_ntc_dt),
        "fold_ntc": partial(fold_ntc, mean_ntc_dt=mean_ntc_dt),
        "contamination_pass": lambda r: r["per_mean_cov_ntc"] <= 1,
        "min_cov_pass": lambda r: r["mean_cov"] >= 50,
        "qc_pass": lambda r: r["contamination_pass"]
        and r["min_cov_pass"]
        and not r["is_negative"],
    }

    for summary_name, summary_func in summary_funcs.items():
        merged_df.insert(
            merged_df.shape[1],
            summary_name,
            merged_df.apply(summary_func, axis=1),  # needs to go row-wise
        )

    # Need to improve this
    expt_summary = (merged_df
                    .query("not is_negative")
                    .groupby("name")
                    .agg(
                        mean_per_mean_cov_ntc=pd.NamedAgg("per_mean_cov_ntc", np.mean),
                        median_per_mean_cov_ntc=pd.NamedAgg("per_mean_cov_ntc", np.median),
                        mean_fold_ntc=pd.NamedAgg("fold_ntc", np.mean),
                        median_fold_ntc=pd.NamedAgg("fold_ntc", np.median),
                        frac_contamination_pass=pd.NamedAgg("contamination_pass", np.mean),
                        frac_min_cov_pass=pd.NamedAgg("min_cov_pass", np.mean),
                        frac_qc_pass=pd.NamedAgg("qc_pass", np.mean)
                    )
                    .reset_index()
                )

    print("Experiment mean Overview:")
    print(expt_summary)
    print("Done.")

    merged_df.to_csv(f"{output_dir}/table.contamiation_qc.csv", index=False)
    ntc_summary_df.to_csv(f"{output_dir}/table.ntc_summaries.csv", index=False)
    expt_summary.to_csv(f"{output_dir}/table.expt_summaries.mean.csv", index=False)

    print_footer(t0)
