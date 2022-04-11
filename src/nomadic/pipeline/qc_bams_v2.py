import pysam
import pandas as pd
pd.options.mode.chained_assignment = None
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
# Load information about alignments from a .bam file
#
# ================================================================


def calc_percent_gc(seq):
    """Calculate GC percentage of a sequence"""

    if not seq:
        return None

    N = len(seq)
    n_gc = len([nt for nt in seq if nt in ["G", "C"]])

    return 100 * n_gc / N


def load_alignment_information(input_bam: str) -> pd.DataFrame:
    """
    Load information about every alignment from an input bam
    file `input_bam`

    """

    @dataclass
    class AlignmentSummary:
        """
        Simple summary of a single alignment
        from a .bam file

        """

        read_id: str
        mapq: int
        flag: int
        query_length: int
        query_alignment_length: int
        mean_qscore: float
        per_gc: float

        @classmethod
        def from_pysam_aligned_segment(cls, pysam_segment: pysam.AlignedSegment):
            """
            Instantiate from a pysam aligned segment

            """
            # Compute mean quality score
            qscores = np.array(pysam_segment.query_qualities)
            mean_qscore = qscores.mean() if qscores.shape else None

            return cls(
                read_id=pysam_segment.query_name,
                mapq=pysam_segment.mapq,
                flag=pysam_segment.flag,
                query_length=pysam_segment.query_length,
                query_alignment_length=pysam_segment.query_alignment_length,
                mean_qscore=mean_qscore,
                per_gc=calc_percent_gc(pysam_segment.query),
            )

    # Iterate over aligned segments in .bam, store results
    with pysam.AlignmentFile(input_bam, "r") as bam:
        results = [
            AlignmentSummary.from_pysam_aligned_segment(segment) for segment in bam
        ]

    return pd.DataFrame(results)


# ================================================================
# Classify reads based on their alignment state
#
# ================================================================


def get_read_mapping_state(flags):
    """
    Classify the alignment of a read, based on its flags

    Five possible outputs:
        - Unmapped (flag == 4)
        - Uniquely mapped (flag in [0, 16])
        - Chimera mapped ([2048, 2064] in flag)
        - Supplementary mapped ([256, 272] in flag)
        - Chimera-Supp mapped (both of above)

    params:
        flags: list of ints
            List of mapping flags associated with this read.

    returns:
        _: str
            A mapping state for the read, draw from
            MapState

    """

    # Handle most common cases first, mapped or unmapped
    if len(flags) == 1:
        flag = flags[0]
        if flag in [0, 16]:
            return "uniq_mapped"
        elif flag == 4:
            return "unmapped"
        else:
            raise ValueError(f"Flag {flag} invalid.")

    # Handle chimera and supplementary
    s = []
    if 2048 in flags or 2064 in flags:
        s.append("chim_mapped")
    if 256 in flags or 272 in flags:
        s.append("supp_mapped")
    if len(s) == 1:
        return s[0]
    elif len(s) > 1:
        # NB: ENFORCING CHIMERA DOMINANT HERE FOR CLARITY
        # Arbitary decision
        return "chim_mapped"
    else:
        raise ValueError(f"Flags {flags} invalid.")


def get_read_state_summary(species, mapping_state, detailed=False):
    """
    Get a single string summary of the read mapping state, taking into
    account `species` and `mapping_state` information

    params
        species: str
            Species to which the read is mapped.
        mapping_state: str
            Mapping state for the read, i.e. unmapped, uniquely
            mapped, chimera mapped, supplementary mapped, or both.
        detailed: bool
            Return more detailed summary

    return
        _: str
            Summary of the species and mapping state.

    """

    # Handle unmapped
    if mapping_state == "unmapped":
        return mapping_state

    return f"{species}_mapped" if not detailed else f"{species}_{mapping_state}"


def reduce_to_read_dataframe(alignments_df):
    """
    Reduce an alignments dataframe to a read-level dataframe,
    annotating mapping state

    Pretty painful function

    """

    # Annotate every alignment with a mapping state
    mapping_state_dt = {
        read_id: get_read_mapping_state(adf["flag"].tolist())
        for read_id, adf in alignments_df.groupby("read_id")
    }
    alignments_df.insert(
        8,
        "mapping_state",
        [mapping_state_dt[read_id] for read_id in alignments_df["read_id"]],
    )

    # Reduce to primary alignment for every read
    read_df = alignments_df.query("flag in [0, 16, 4]")

    # Annotate with summary of species + mapping state
    read_df.insert(
        9,
        "primary_state",
        [
            get_read_state_summary(s, m)
            for s, m in zip(read_df["species"], read_df["mapping_state"])
        ],
    )
    read_df.insert(
        10,
        "secondary_state",
        [
            get_read_state_summary(s, m, detailed=True)
            for s, m in zip(read_df["species"], read_df["mapping_state"])
        ],
    )

    # Sort levels
    # NB: For now I am dropping both chimera and supplementary, too hard to view
    primary_levels = ["pf_mapped", "hs_mapped", "unmapped"][::-1]
    read_df.loc[:, "primary_state"] = pd.Categorical(
        values=read_df["primary_state"], categories=primary_levels, ordered=True
    )

    secondary_levels = ["unmapped"] + [
        f"{s}_{m}"
        for s in ["pf", "hs"]
        for m in ["uniq_mapped", "chim_mapped", "supp_mapped"]
    ][::-1]
    read_df.loc[:, "secondary_state"] = pd.Categorical(
        values=read_df["secondary_state"], categories=secondary_levels, ordered=True
    )

    return read_df


# ================================================================
# Plotting functionality
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

    def get_histogram_dataframe(self):
        """
        This requires `.plot_histogram()`
        is run first

        """

        hist_df = pd.DataFrame(self.hist_info[0].transpose())
        hist_df.columns = [n for n, _ in self.grps]
        hist_df.index = [
            f"({low}-{high}]"
            for low, high in zip(self.hist_info[1][:-1], self.hist_info[1][1:])
        ]

        return hist_df

    def plot_histogram(self, title=None, output_path=None):
        """
        Plot from histograms

        """
        fig, ax = plt.subplots(1, 1, figsize=(10, 4))

        # Prepare labels
        @dataclass
        class GroupInfo:
            name: str  # Name of group
            n: int  # Number of reads in group
            mu: float  # Mean of group statistic
            N: int  # Total reads

            @classmethod
            def from_df(cls, name, df, stat, N):
                return cls(name=name, n=df.shape[0], mu=df[stat].mean(), N=N)

            def __repr__(self):
                return f"{self.name} ($n=${self.n}, {100*self.n/self.N:.1f}%)"

        labels = [
            GroupInfo.from_df(name=name, df=gdf, stat=self.stat, N=self.N)
            for name, gdf in self.grps
        ]

        # Plot
        self.hist_info = ax.hist(
            [df[self.stat] for _, df in self.grps],
            label=[str(label) for label in labels],
            color=self.colors,
            bins=self.bins,
            histtype="stepfilled",
            ls="solid",
            ec="black",
            lw=0.5,
            stacked=True,
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

        # Define input bams
        pf_bam_path = f"{bam_dir}/{barcode}.{pf_reference.name}.final.sorted.bam"
        hs_bam_path = f"{bam_dir}/{barcode}.{hs_reference.name}.final.sorted.bam"

        # Load p.f. alignments
        print("Loading data...")
        pf_alignments_df = load_alignment_information(pf_bam_path)
        pf_alignments_df.insert(0, "species", "pf")
        pf_alignments_df.query(
            "flag != 4", inplace=True
        )  # these have been remapped to H.s.

        # Load h.s. alignments
        hs_alignments_df = load_alignment_information(hs_bam_path)
        hs_alignments_df.insert(0, "species", "hs")

        # Combine all alignments
        alignments_df = pd.concat([pf_alignments_df, hs_alignments_df])

        # Produce a read-level data frame
        print("Processing...")
        read_df = reduce_to_read_dataframe(alignments_df)

        # Prepare plotter
        plotter = ReadHistogramPlotter(read_df)

        # Prepare colors
        pf_colors = sns.color_palette("Blues_r", 3)
        hs_colors = sns.color_palette("Reds_r", 3)
        unmapped_color = (0.75, 0.75, 0.75)

        primary_colors = [unmapped_color, hs_colors[0], pf_colors[0]]
        secondary_colors = pf_colors + hs_colors + [unmapped_color]

        color_sets = {
            "primary_state": primary_colors,
            "secondary_state": secondary_colors[::-1],
        }

        # Iterate over statistics, states, and plot
        print("Plotting...")
        for state, colors in color_sets.items():
            plotter.set_groups(state, colors)
            size_df = plotter.get_group_size_dataframe()
            size_df.to_csv(f"{output_dir}/table.size.{state}.csv", index=False)

            for histogram_stat in histogram_stats:
                plotter.set_histogram_stats(**histogram_stat.__dict__)

                # Plot
                plotter.plot_histogram(
                    title=barcode,
                    output_path=f"{output_dir}/plot.{histogram_stat.stat}.{state}.png",
                )

                # Write to data frame
                hist_df = plotter.get_histogram_dataframe()
                hist_df.to_csv(
                    f"{output_dir}/table.bin_counts.{histogram_stat.stat}.{state}.csv",
                    index=True,
                )
        print(f"Output directory: {output_dir}")
        print("Done.")
        print("")

    print_footer(t0)
