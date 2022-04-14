import pandas as pd
import numpy as np
import random

from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.parsing import build_parameter_dict

import matplotlib.pyplot as plt
import seaborn as sns


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


class BalancePlotter:

    scale_size = 0.25

    def __init__(self, sample_ids, gene_names, values):
        """
        Plot a statistics `values` across `sample_ids` and `genes`

        """

        # Store instance variables
        self.sample_ids = sample_ids
        self.gene_names = gene_names
        self.values = values

        # Create data frame
        self.df = pd.DataFrame(
            {"sample_id": sample_ids, "gene_name": gene_names, "value": values}
        )

        # UNIQUE
        # sample_ids
        self.unique_sample_ids = np.unique(self.sample_ids)
        self.n_sample_ids = self.unique_sample_ids.shape[0]
        # gene_names
        self.unique_gene_names = np.unique(self.gene_names)
        self.n_gene_names = self.unique_gene_names.shape[0]

        # Create sample_id axis map
        self._sample_axis_map()

    def _sample_axis_map(self):
        """
        Create a map from sample id to axis position

        """

        # Prepare ticks and limits
        self.sample_pos = np.arange(self.n_sample_ids)
        self.sample_pos_min = self.sample_pos.min()
        self.sample_pos_max = self.sample_pos.max()

        # Map
        self.id_to_pos = {
            sample_id: pos
            for sample_id, pos in zip(self.unique_sample_ids, self.sample_pos)
        }
        self.pos_to_id = {  # safer to reverse map, in case needed
            pos: sample_id for sample_id, pos in self.id_to_pos.items()
        }

    def get_sample_axis_pos(self, sample_id, jitter=True):
        """
        Return the axis position of a sample

        """

        y = self.id_to_pos[sample_id]
        if jitter:
            y += (random.random() - 0.5) * self.scale_size * 1.8

        return y

    def set_gene_color_pal(self, pal):
        """
        Set color palette for the genes

        """
        self.gene_col = dict(
            zip(self.unique_gene_names, sns.color_palette(pal, self.n_gene_names))
        )

    def plot(self, output_path=None):
        """
        Create the balance plot

        """

        # Prepare canvas
        fig, ax = plt.subplots(figsize=(5, self.scale_size * self.n_sample_ids))

        # Plot
        for gene, gene_df in self.df.groupby("gene_name"):
            ax.scatter(
                y=gene_df["sample_id"].apply(self.get_sample_axis_pos),
                x=gene_df["value"],
                color=self.gene_col[gene],
                label=gene,
                ec="black",
                alpha=0.8,
            )

        # Scale
        ax.set_xscale("log")

        # Sample limits, ticks and labels
        ax.set_ylim(self.sample_pos_min - 1, self.sample_pos_max + 1)
        ax.set_yticks(self.sample_pos)
        ax.set_yticklabels(map(lambda p: self.pos_to_id[p], self.sample_pos))

        # Value label
        ax.set_xlabel("No. Reads\nOverlapping Target")

        # Invert yaxis
        ax.invert_yaxis()

        # Ticks
        ax.grid(ls="dotted", axis="x", zorder=-10)
        ax.axvline(x=100, color="black", zorder=-10)

        # Delineate barcodes
        improve_delin = True
        if improve_delin:
            for j in range(self.n_sample_ids)[::2]:
                ax.axhline(
                    j, lw=self.scale_size * 50, color="lightgrey", alpha=0.5, zorder=-10
                )

        # Legend
        ax.legend(loc="lower left", bbox_to_anchor=(0, 1), ncol=3)

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

        # Reduce to plotting information
        plot_df = merged_df.query("overlap == @overlap_type")

        # Pivot for heatmap
        pivot_df = plot_df.pivot(
            index="sample_id", columns="gene_name", values="reads_mapped"
        )

        # Plot heatmap
        plot_dataframe_heat(
            pivot_df,
            cmap="Wistia",
            output_path=f"{output_dir}/plot.reads_mapped.heatmap.{overlap_type}.pdf",
        )

        # Plot stripplot
        plotter = BalancePlotter(
            sample_ids=plot_df["sample_id"],
            gene_names=plot_df["gene_name"],
            values=plot_df["reads_mapped"],
        )
        plotter.set_gene_color_pal("Spectral")
        plotter.plot(
            output_path=f"{output_dir}/plot.reads_mapped.stripplot.{overlap_type}.pdf"
        )

    print_footer(t0)
