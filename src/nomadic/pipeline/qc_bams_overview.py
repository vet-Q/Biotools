import pandas as pd
import numpy as np
from dataclasses import dataclass
from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.pipeline.qc_bams import MappingStatesAndColors, histogram_stats

import matplotlib.pyplot as plt
import seaborn as sns


# ================================================================
# Load dataframes across all barcodes in the experiment
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


class JointHistogramPlotter:
    def __init__(self, df, columns, colors):
        """
        Plot a histogram  across all experiments

        Perhaps could inheret from `ReadHistogramPlotter`,
        not so critical though

        """
        self.df = df
        self.N = self.df[columns].sum().sum()
        self.columns = columns
        self.colors = colors

    def set_histogram_stats(
        self, stat, name, min_val, max_val, intv, major_loc=None, minor_loc=None
    ):
        """
        Set histogram statistics

        """
        self.stat = stat
        self.name = name
        self.min_val = min_val
        self.max_val = max_val
        self.intv = intv
        self.bins = np.arange(min_val, max_val + intv, intv)
        self.major_loc = major_loc
        self.minor_loc = minor_loc

    def plot_histogram(self, title=None, output_path=None):
        """
        Plot histogram

        """

        # Prepare labels
        @dataclass
        class GroupInfo:
            name: str  # Name of group
            n: int  # Number of reads in group
            # mu: float  # Mean of group statistic
            N: int  # Total reads

            @classmethod
            def from_column(cls, name, df, column, N):

                return cls(name=name, n=df[column].sum(), N=N)

            def __repr__(self):
                return f"{self.name} ($n=${self.n}, {100*self.n/self.N:.1f}%)"

        labels = [
            GroupInfo.from_column(name=c, df=self.df, column=c, N=self.N)
            for c in self.columns
        ]

        fig, ax = plt.subplots(1, 1, figsize=(10, 4))

        # Plot
        ax.stackplot(
            self.df["bin_lower"],
            self.df[self.columns].transpose(),
            labels=labels,
            step="mid",
            colors=self.colors,
            ec="black",
            lw=0.8,
        )

        # Labels
        ax.set_xlabel(self.name)
        ax.set_ylabel("No. reads")
        if title is not None:
            ax.set_title(title, loc="left")

        # Limits & ticks
        ax.set_xlim(0 - self.intv, self.max_val + self.intv)
        if self.minor_loc:
            ax.xaxis.set_minor_locator(plt.MultipleLocator(self.minor_loc))
        if self.major_loc:
            ax.xaxis.set_major_locator(plt.MultipleLocator(self.major_loc))

        # Legend
        ax.legend(title="Reads", loc="upper right")

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

        # Tall to wide
        joint_df = (joint_df
                    .pivot(index="barcode", columns="group", values="n_reads")
                    .reset_index()
                )

        # Merge metadata
        merged_df = pd.merge(left=joint_df, right=params["metadata"], on="barcode")
        assert (
            joint_df.shape[0] == merged_df.shape[0]
        ), "Barcodes lost during merge with metadata."

        # Write
        merged_df.to_csv(f"{output_dir}/table.mapping.{state}.csv", index=False)

        # Isolate to key columns for plotting
        plot_df = merged_df[msc.level_sets[state][::-1]]
        plot_df.index = merged_df["sample_id"]

        # Barplot
        barplot_states(
            plot_df,
            colors=msc.color_sets[state][::-1],
            show_per="pf_uniq_mapped" if state == "secondary_state" else "pf_mapped",
            output_path=f"{output_dir}/plot.mapping.{state}.pdf",
        )

        # Now plot histograms
        for histogram_stat in histogram_stats:
            # Load data
            joint_df = combine_barcode_dataframes(
                csv_filename=f"table.bin_counts.{histogram_stat.stat}.{state}.csv",
                script_dir=script_dir,
                params=params,
            )

            # Aggregate
            sum_df = joint_df.groupby(["bin_lower", "bin_higher"]).sum().reset_index()

            # Plot histogram
            plotter = JointHistogramPlotter(
                sum_df, msc.level_sets[state], msc.color_sets[state]
            )
            plotter.set_histogram_stats(**histogram_stat.__dict__)
            plotter.plot_histogram(
                output_path=f"{output_dir}/plot.hist.{histogram_stat.stat}.{state}.pdf"
            )
    print_footer(t0)
