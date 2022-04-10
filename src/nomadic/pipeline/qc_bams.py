import pysam
import pandas as pd
import numpy as np
from dataclasses import dataclass
import matplotlib.pyplot as plt
import seaborn as sns

from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.references import (
    PlasmodiumFalciparum3D7,
    HomoSapiens,
)


# ================================================================
# Extract aligment (/read) information from a .bam file
#
# ================================================================


def calc_percent_gc(seq):
    """Calculate GC percentage of a sequece"""

    if not seq:
        return None

    N = len(seq)
    n_gc = len([nt for nt in seq if nt in ["G", "C"]])

    return 100 * n_gc / N


def extract_read_information(input_bam):
    """Extract a dataframe containing information about reads"""

    # Define a small class to hold information about reads from .bam
    @dataclass
    class ReadInfo:
        mapq: int
        flag: int
        query_length: int
        mean_qscore: float
        per_gc: float

    # Iterate over reads in .bam, store results
    results = []
    with pysam.AlignmentFile(input_bam, "r") as bam:
        for alignment in bam:
            # Compute mean quality score
            qscores = np.array(alignment.query_qualities)
            mean_qscore = qscores.mean() if qscores.shape else None

            # Store
            read_info = ReadInfo(
                mapq=alignment.mapq,
                flag=alignment.flag,
                query_length=alignment.query_length,
                mean_qscore=mean_qscore,
                per_gc=calc_percent_gc(alignment.query),
            )
            results.append(read_info)

    return pd.DataFrame(results)


# ================================================================
# Read grouping
#
# ================================================================


class ReadGrouping:
    def __init__(self, df, queries):
        """
        Create a set of groupings for a dataframe `df` based on
        a set of query strings in `queries`

        Each query is applied as `df.query(query)`, see
        `._create_groupings()`

        The result is to split the input dataframe, `df`,
        into a set of dictionaries `dfs`, where the key
        is the name of the query, and the value is the
        dataframe resulting from that query.

        """
        self.queries = queries
        self.n_grps = len(self.queries)
        self.df = df
        self.dfs = self._create_groupings()

    def _create_groupings(self):
        """
        Apply a set of `queries` to a dataframe `df`,
        returning a dictionary of queries

        """
        dfs = {l: self.df.query(q) for l, q in self.queries.items()}

        return dfs

    def get_group_sizes(self):
        """
        Get the sizes of the groups, returned as a dataframe

        """
        return pd.DataFrame(
            [(l, df.shape[0]) for l, df in self.dfs.items()],
            columns=["read_group", "n_reads"],
        )


# ================================================================
# Plotting
#
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
histogram_stats = [read_length, gc, qscore]


class ReadGroupingHistogram:
    def __init__(self, grouping):
        self.grouping = grouping
        self.stat = None
        self.name = None
        self.min_val = None
        self.max_val = None
        self.intv = None
        self.bins = None
        self.minor_loc = None
        self.major_loc = None
        self.colors = None

    def set_plot_statistic(
        self,
        stat,
        name,
        min_val,
        max_val,
        intv,
        minor_loc=None,
        major_loc=None,
        colors=None,
    ):
        """Set statistic and other plotting parameters"""

        self.stat = stat
        assert (
            self.stat in self.grouping.df.columns
        ), f"Statistic {self.stat} not in grouping columns: {', '.join(self.grouping.df.columns)}."
        self.name = name
        self.min_val = min_val
        self.max_val = max_val
        self.intv = intv
        self.bins = np.arange(min_val, max_val + intv, intv)
        self.minor_loc = minor_loc
        self.major_loc = major_loc
        if colors is not None:
            self.colors = colors
        else:
            self.colors = sns.color_palette("Set1", self.grouping.n_grps)

        return None

    def calc_statistic_bin_counts(self):
        """Calculate histogram bin counts for a given statistic"""

        assert self.stat, "`.set_plot_statistic()` must be run first."

        dt = {"stat": [self.stat] * len(self.bins[1:]), "bins": self.bins[1:]}
        for l, df in self.grouping.dfs.items():
            counts, _ = np.histogram(df[self.stat], bins=self.bins)
            dt[l] = counts

        return pd.DataFrame(dt)

    def plot(self, title=None, annot_med=False, output_path=None):
        """Plot histogram of target statistics over read groupings"""

        fig, ax = plt.subplots(1, 1, figsize=(10, 4))

        # Prepare for plotting
        stats = {
            l: (df.shape[0], df[self.stat].median())
            for l, df in self.grouping.dfs.items()
        }
        N = sum([s[0] for l, s in stats.items()])
        color_dt = {l: c for l, c in zip(self.grouping.queries, self.colors)}

        # Plot
        ax.hist(
            [df[self.stat] for l, df in self.grouping.dfs.items()],
            label=[f"{l} ($n=${s[0]}, {100*s[0]/N:.1f}%)" for l, s in stats.items()],
            color=[c for l, c in color_dt.items()],
            ls="solid",
            histtype="step",
            bins=self.bins,
            stacked=False,
        )

        # Limits & Ticks
        ax.set_xlim(0 - self.intv, self.max_val + self.intv)
        if self.minor_loc:
            ax.xaxis.set_minor_locator(plt.MultipleLocator(self.minor_loc))
        if self.major_loc:
            ax.xaxis.set_major_locator(plt.MultipleLocator(self.major_loc))

        # Summary statistics
        for l, (n, m) in stats.items():
            ax.axvline(x=m, color=color_dt[l], lw=1, ls="dashed")
            if annot_med:
                ax.annotate(
                    xy=(m, ax.get_ylim()[1]),
                    ha="right",
                    va="top",
                    rotation=90,
                    text="$\~{\mu}=$%.0f " % m,
                    color=color_dt[l],
                )

        # Legend
        ax.set_xlabel(self.name)
        ax.set_ylabel("Count")
        ax.legend(title="Reads", loc="upper right")
        ax.set_title(title, loc="left")

        # Save
        if output_path is not None:
            fig.savefig(output_path, bbox_inches="tight", pad_inches=0.5, dpi=300)
            plt.close(fig)


# ================================================================
# Main script, run from `cli.py`
#
# ================================================================


def main(expt_dir, config, barcode):
    """
    Create a series of histogram summaries of reads
    within .bam files

    """

    # PARSE INPUTS
    script_descrip = "NOMADIC: Quality control analysis of .bam files"
    t0 = print_header(script_descrip)
    script_dir = "qc-bams"
    params = build_parameter_dict(expt_dir, config, barcode)

    # Focus on a single barcode, if specified
    if "focus_barcode" in params:
        params["barcodes"] = [params["focus_barcode"]]

    # Define reference genomes
    pf_reference = PlasmodiumFalciparum3D7()
    hs_reference = HomoSapiens()

    # Define queries
    queries = {
        "mapping": {
            "Unmapped": "species == 'Hs' and flag == 4",
            "Mapped H.s.": "species == 'Hs' and flag != 4",
            "Mapped P.f.": "species == 'Pf'",
        },
        "flag": {
            flag: f"species == 'Pf' and flag == {flag}" for flag in [2064, 2048, 16, 0]
        },
    }
    queries_colors = {
        "mapping": ["#a9a9a9", "#e41a1c", "#377eb8"],
        "flag": sns.color_palette("Paired", 8)[:2]
        + sns.color_palette("Paired", 8)[-2:],
    }

    # ITERATE
    print("Iterating over barcodes...")
    for barcode in params["barcodes"]:
        print("." * 80)
        print(f"Barcode: {barcode}")
        print("." * 80)

        # Define directories
        barcode_dir = f"{params['barcodes_dir']}/{barcode}"
        bam_dir = f"{barcode_dir}/bams"
        output_dir = produce_dir(barcode_dir, script_dir)
        print(barcode_dir)
        print(output_dir)

        # Define input bams
        pf_bam_path = f"{bam_dir}/{barcode}.{pf_reference.name}.final.sorted.bam"
        hs_bam_path = f"{bam_dir}/{barcode}.{hs_reference.name}.final.sorted.bam"

        # Load Pf reads
        print("Loading P.f. reads...")
        pf_read_df = extract_read_information(pf_bam_path)
        pf_read_df.insert(0, "species", "Pf")
        pf_read_df.query(
            "flag != 4", inplace=True
        )  # drop unmapped, because have been remapped to H.s.

        # Load Hs reads
        print("Loading H.S. reads...")
        hs_read_df = extract_read_information(hs_bam_path)
        hs_read_df.insert(0, "species", "Hs")

        # Combine
        read_df = pd.concat([pf_read_df, hs_read_df])

        # Iterate over groupings
        print("Plotting...")
        for query_name, query_groups in queries.items():

            # Create read grouping and plotter
            grouping = ReadGrouping(read_df, query_groups)
            group_df = grouping.get_group_sizes()
            group_df.to_csv(f"{output_dir}/table.sizes.{query_name}.csv", index=False)

            # Instantiate a histogram plotter for the group
            plotter = ReadGroupingHistogram(grouping)

            # Create and save plots
            for s in histogram_stats:
                plotter.set_plot_statistic(
                    **s.__dict__, colors=queries_colors[query_name]
                )
                plotter.plot(
                    annot_med=True,
                    output_path=f"{output_dir}/plot.{s.stat}.{query_name}.png",
                )
                hist_df = plotter.calc_statistic_bin_counts()
                hist_df.to_csv(
                    f"{output_dir}/table.bin_counts.{s.stat}.{query_name}.csv",
                    index=False,
                )
        print("Done.")
        print("")

    print_footer(t0)
