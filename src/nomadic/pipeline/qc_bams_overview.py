import os
import pandas as pd
from nomadic.lib.parsing import build_parameter_dict
from nomadic.pipeline.qc_bams_v2 import MappingStatesAndColors
from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.parsing import build_parameter_dict

import matplotlib.pyplot as plt
import seaborn as sns


# ================================================================
# Load dataframes across all barcodes in the experiment
#
# ================================================================


def combine_barcode_dataframes(script_dir, csv_filename, params):
    """
    Combine a particular dataframe across all barcodes

    """

    dfs = []
    for barcode in params["barcodes"]:
        csv_path = f"{params['barcodes_dir']}/{barcode}/{script_dir}/{csv_filename}"
        df = pd.read_csv(csv_path)
        df.insert(0, "barcode", barcode)
        dfs.append(df)

    return pd.concat(dfs)


# ================================================================
# Plotting
#
# ================================================================


def barplot_states(
    df,
    colors=None,
    size_scale=0.25,
    pwr=3,
    show_per=None,
    improve_delin=True,
    output_path=None,
):
    """
    Make a barplot of the mapping states in `df`

    """

    n_barcodes = df.shape[0]
    n_states = df.shape[1]
    if colors is None:
        colors = sns.color_palette("viridis", n_states)

    # Prepare canvas
    fig, ax = plt.subplots(figsize=(5, size_scale * n_barcodes))

    # Plot
    df.plot(kind="barh", stacked=True, ec="black", width=0.8, color=colors, ax=ax)

    # Axis inversion
    ax.invert_yaxis()

    # Grid
    ax.set_axisbelow(True)
    ax.grid(axis="x", ls="dotted", alpha=0.5)

    # Ticks
    ax.xaxis.set_major_formatter(
        plt.FuncFormatter(lambda val, _: f"{int(val/(10**pwr))}")
    )

    # Labels
    ax.set_xlabel(f"No. Reads [$x10^{pwr}]$")
    ax.set_ylabel("")

    # Legend
    ax.legend(loc="lower left", bbox_to_anchor=(-0.05, 1), ncol=3)

    # Delineate barcodes
    if improve_delin:
        for j in range(n_barcodes)[::2]:
            ax.axhline(j, lw=10, color="lightgrey", alpha=0.5, zorder=-10)

    if show_per is not None:
        ix = list(df.columns).index(show_per)
        show_per_color = colors[ix]
        for i, (_, row) in enumerate(df.iterrows()):
            ax.annotate(
                xy=(ax.get_xlim()[1], ax.get_yticks()[i]),
                xycoords="data",
                ha="left",
                va="center",
                color=show_per_color,
                text=f"  {100*row[show_per]/row.sum():.1f}%",
            )

    # Save
    if output_path is not None:
        fig.savefig(output_path, bbox_inches="tight", pad_inches=0.5, dpi=300)
        plt.close(fig)


# ================================================================
# Main script, run from `cli.py`
#
# ================================================================


def main(expt_dir, config):
    """
    Produce a summary of .bam QC across all barcodes within an expeirment

    """
    # PARSE INPUTS
    script_descrip = "NOMADIC: Overview of .bam file QC across all barcodes"
    t0 = print_header(script_descrip)
    script_dir = "qc-bams"
    params = build_parameter_dict(expt_dir, config, barcode=None)
    output_dir = produce_dir(params["nomadic_dir"], script_dir)

    # Instantiate mapping states and colors
    msc = MappingStatesAndColors()

    # Iterate over states and produce barplots
    states = ["primary_state", "secondary_state"]
    for state in states:
        # Load data
        joint_df = combine_barcode_dataframes(
            csv_filename=f"table.size.{state}.csv", script_dir=script_dir, params=params
        )

        # Merge metadata
        # - Will need a standardised ID for every sample in the metadata for this to work
        merged_df = pd.merge(left=joint_df, right=params["metadata"], on="barcode")
        assert (
            joint_df.shape[0] == merged_df.shape[0]
        ), "Barcodes lost during merge with metadata."

        # Write
        merged_df.to_csv(f"{output_dir}/table.mapping.{state}.csv")

        # Wide for plotting
        wide_df = merged_df.pivot_table(index="name", columns="group", values="n_reads")
        wide_df = wide_df[msc.level_sets[state][::-1]]

        # Plot
        barplot_states(
            wide_df,
            colors=msc.color_sets[state][::-1],
            show_per="pf_uniq_mapped" if state == "secondary_state" else "pf_mapped",
            output_path=f"{output_dir}/plot.mapping.{state}.pdf",
        )

    print_footer(t0)
