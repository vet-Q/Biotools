
import pandas as pd
import numpy as np

from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.parsing import build_parameter_dict

import matplotlib.pyplot as plt


# ================================================================
# Load data
#
# ================================================================


def combine_barcode_dataframes(csv_filename, script_dir, params):
    """
    Combine a particular dataframe across all barcodes

    """

    dfs = []
    for barcode in params["barcodes"]:
        csv_path = f"{params['barcodes_dir']}/{barcode}/{script_dir}/{csv_filename}"
        df = pd.read_csv(csv_path)
        # df.insert(0, "barcode", barcode)
        dfs.append(df)

    return pd.concat(dfs)


# ================================================================
# Plotting
#
# ================================================================


def plot_dataframe_heat(
    df,
    title=None,
    cmap="bwr",
    cbar=None,
    as_log=False,
    label_values=True,
    output_path=None,
    **kwargs,
):
    """
    Create a `plt.imshow` like visualisation of a data frame

    """

    # Set figure size proportional to data frame
    r, c = df.shape
    rescale = 0.6
    fig, ax = plt.subplots(1, 1, figsize=(c * rescale, r * rescale))

    # Plot
    cax = ax.matshow(np.log(df) if as_log else df, cmap=cmap, **kwargs)

    # Colorbar
    if cbar is not None:
        fig.colorbar(
            cax,
            shrink=0.85 if c > 10 else 0.8,
            pad=0.01 if c > 10 else 0.05,
            label=cbar,
        )

    # Axis, Ticks, Grid
    ax.xaxis.set_label_position("top")
    ax.set_yticks(range(r))
    ax.set_xticks(range(c))
    ax.set_yticks(np.arange(0.5, r + 0.5), minor=True)
    ax.set_xticks(np.arange(0.5, c + 0.5), minor=True)
    ax.set_yticklabels(df.index)
    ax.set_xticklabels(df.columns, rotation=90)
    ax.grid(which="minor", ls="dotted")

    # Optional title
    if title is not None:
        ax.set_title(title)

    # Text labels
    if label_values:
        tau = np.quantile(a=df, q=0.2)
        for i in range(r):
            for j in range(c):
                v = df.iloc[i, j]
                if 0 <= v <= 1:
                    t = "%.02f" % v
                else:
                    t = "%.0f" % v
                if cmap in ["bwr"]:
                    ax.text(
                        j,
                        i,
                        t,
                        va="center",
                        ha="center",
                        size=7,
                        color="black" if v >= tau else "white",
                    )
                else:
                    ax.text(j, i, t, va="center", ha="center", size=7, color="black")

    # Optionally save figure
    if output_path is not None:
        fig.savefig(output_path, bbox_inches="tight", pad_inches=0.5, dpi=300)
        plt.close(fig)

    return None


# ================================================================
# Main script, run from `cli.py`
#
# ================================================================


def main(expt_dir, config):
    """
    Overview of coverage across amplicons

    """

    # PARSE INPUTS
    script_descrip = "NOMADIC: Overview of target extraction."
    t0 = print_header(script_descrip)
    script_dir = "target-extraction"
    params = build_parameter_dict(expt_dir, config, barcode=None)
    output_dir = produce_dir(params["nomadic_dir"], script_dir)

    # LOAD DATA
    joint_df = combine_barcode_dataframes(
        csv_filename="table.extraction.summary.csv",
        script_dir="target-extraction",
        params=params,
    )

    # Merge with metadata
    merged_df = pd.merge(left=joint_df, right=params["metadata"], on="barcode")

    # Write
    merged_df.to_csv(f"{output_dir}/table.target_coverage.overview.csv", index=False)

    # Iterate over overlap types and plot
    overlap_types = ["complete", "any"]
    for overlap_type in overlap_types:
        # Pivot
        pivot_df = merged_df.query("overlap == @overlap_type").pivot(
            index="sample_id", columns="gene_name", values="reads_mapped"
        )

        # Plot
        plot_dataframe_heat(
            pivot_df,
            cmap="Wistia",
            output_path=f"{output_dir}/plot.target_coverage.{overlap_type}.pdf",
        )

    print_footer(t0)
