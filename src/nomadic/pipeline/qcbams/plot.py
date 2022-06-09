import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from functools import partial
from dataclasses import dataclass

# ================================================================
# Map alignment states to colors
#
# ================================================================


class MappingStatesAndColors:
    """
    Class to hold mapping states and their associated colors for
    plotting

    This class is strongly coupled to `get_read_mapping_state()`
    and `get_read_state_summary()`

    """

    # Species / states
    species = ["pf", "hs"]
    mapping_states = ["uniq_mapped", "chim_mapped", "supp_mapped"]
    unmapped_state = ["unmapped"]

    def __init__(self):
        self._prepare_levels()
        self._prepare_colors()

    def _prepare_levels(self):
        self.primary_levels = (
            self.unmapped_state + [f"{s}_mapped" for s in self.species][::-1]
        )
        self.secondary_levels = (
            self.unmapped_state
            + [f"{s}_{m}" for s in self.species for m in self.mapping_states][::-1]
        )
        self.level_sets = {
            "primary_state": self.primary_levels,
            "secondary_state": self.secondary_levels,
        }

    def _prepare_colors(self):
        # Colors
        self.pf_colors = sns.color_palette("Blues_r", len(self.mapping_states))
        self.hs_colors = sns.color_palette("Reds_r", len(self.mapping_states))
        self.unmapped_color = [(0.75, 0.75, 0.75)]

        self.primary_colors = self.unmapped_color + [
            self.hs_colors[0],
            self.pf_colors[0],
        ]
        self.secondary_colors = (
            self.unmapped_color + self.hs_colors[::-1] + self.pf_colors[::-1]
        )
        self.color_sets = {
            "primary_state": self.primary_colors,
            "secondary_state": self.secondary_colors,
        }


# ================================================================
# Store metadata needed for plotting histograms across different
# statistics
# ================================================================


@dataclass
class HistogramStatistic:
    """Information necessary to plot a histogram of a given statistic"""

    stat: str
    name: str
    min_val: int
    max_val: int
    intv: int
    minor_loc: int


read_length = HistogramStatistic(
    stat="query_length",
    name="Read length (bp)",
    min_val=0,
    max_val=6000,
    intv=50,
    minor_loc=500,
)
gc = HistogramStatistic(
    stat="per_gc",
    name="GC (%)",
    min_val=0,
    max_val=100,
    intv=1,
    minor_loc=10,
)
qscore = HistogramStatistic(
    stat="mean_qscore",
    name="Mean Q-score",
    min_val=0,
    max_val=60,
    intv=0.5,
    minor_loc=5,
)
HISTOGRAM_STATS = [read_length, gc, qscore]


# ================================================================
# Plot histograms
#
# ================================================================


class ReadHistogramPlotter:
    def __init__(self, read_df):
        """
        Plot a histogram for reads found in `read_df`

        """
        self.read_df = read_df
        self.N = read_df.shape[0]

    def set_groups(self, group_column, colors):
        """
        Group the reads based on a single column

        """

        self.grps = self.read_df.groupby(group_column, sort=False)
        self.n_grps = len(self.grps)
        self.colors = colors
        self.columns = [n for n, _ in self.grps]

    def get_group_size_dataframe(self):
        """
        Get sizes of the different groups,
        requires groups are set

        """
        size_df = pd.DataFrame(self.grps.size()).reset_index()
        size_df.columns = ["group", "n_reads"]

        return size_df

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

    def create_histogram_dataframe(self):
        """
        Create a data frame containing information
        necessary to plot histogram


        """

        # Define histogram function
        hist_func = partial(np.histogram, bins=self.bins)

        # Call on groups of reads
        self.hist_df = pd.DataFrame(
            {grp_name: hist_func(df[self.stat])[0] for grp_name, df in self.grps}
        )

        # Insert bin boundaries
        self.hist_df.insert(0, "bin_lower", self.bins[:-1])
        self.hist_df.insert(1, "bin_higher", self.bins[1:])

    def write_histogram_dataframe(self, output_path):
        """
        Write the histogram dataframe to `output_path`

        """
        self.hist_df.to_csv(output_path, index=False)

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
            GroupInfo.from_column(name=c, df=self.hist_df, column=c, N=self.N)
            for c in self.columns
        ]

        fig, ax = plt.subplots(1, 1, figsize=(10, 4))

        # Plot
        ax.stackplot(
            self.hist_df["bin_lower"],
            self.hist_df[self.columns].transpose(),
            labels=labels,
            step="mid",
            colors=self.colors,
            ec="black",
            lw=0.5,
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
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], title="Reads", loc="upper right")

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
            lw=0.5,
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
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], title="Reads", loc="upper right")

        # Save
        if output_path is not None:
            fig.savefig(output_path, bbox_inches="tight", pad_inches=0.5, dpi=300)
            plt.close(fig)


# ================================================================
# Create barplots
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
    def get_spacing(spacings, max_val, max_ticks=8):
        for spacing in spacings:
            ticks = np.arange(0, max_val + spacing, spacing)
            if len(ticks) < max_ticks:
                break
        return spacing
    spacing = get_spacing(
        spacings=np.array([1, 2, 5, 10, 20, 50])*10**pwr,
        max_val=df.sum(1).max()
    )
    ax.xaxis.set_major_locator(plt.MultipleLocator(spacing))
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda val, _: f"{int(val/(10**pwr))}"))

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

