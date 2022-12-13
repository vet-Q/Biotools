import click
import pysam
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from functools import partial
from collections import namedtuple
from dataclasses import dataclass
from matplotlib.gridspec import GridSpec

from nomadic.lib.process_gffs import load_gff, add_gff_fields
from nomadic.lib.references import reference_collection
from nomadic.truthset.fasta import load_haplotype_from_fasta
from nomadic.lib.statistics import (
    calc_sliding_percentGC,
    get_homopolymer_runs,
    get_array_encoding,
)


# --------------------------------------------------------------------------------
# Aesthetic parameters
#
# --------------------------------------------------------------------------------







# --------------------------------------------------------------------------------
# Data loading
#
# --------------------------------------------------------------------------------


@dataclass
class AlignmentPlotData:
    read_id: str
    ref_start: int
    ref_end: int
    ref_name: str
    query_alignment_length: int
    mean_qscore: float
    forward: bool

    @classmethod
    def from_pysam_aligned_segment(cls, pysam_segment: pysam.AlignedSegment):
        """
        Extract alignment data for plotting from an aligned `pysam` segment

        """
        # Compute mean quality score
        qscores = np.array(pysam_segment.query_qualities)
        mean_qscore = qscores.mean() if qscores.shape else None

        return cls(
            read_id=pysam_segment.query_name,
            ref_start=pysam_segment.reference_start,
            ref_end=pysam_segment.reference_end,
            ref_name=pysam_segment.reference_name,
            query_alignment_length=pysam_segment.query_alignment_length,
            mean_qscore=mean_qscore,
            forward=pysam_segment.is_forward,
        )


# --------------------------------------------------------------------------------
# Individual plotting classes
#
# --------------------------------------------------------------------------------


class GffPlotter:

    gff_features = ["protein_coding_gene", "CDS"]

    def __init__(self, gff, chrom, start, end, amplicon_region=None):
        """
        Plot information from a GFF over a defined
        region

        """
        self.gff = gff
        self.chrom = chrom
        self.start = start
        self.end = end
        self.amplicon_region=amplicon_region

        self.plot_gff = self._get_relevant_gff()

    def _get_relevant_gff(self):
        """Get the relevant proportion of the .gff for plotting"""

        qry = f"(seqname == '{self.chrom}')"
        qry += f" and ({self.start} <= start <= {self.end}"
        qry += f" or {self.start} <= end <= {self.end}"
        qry += f" or (start <= {self.start} and {self.end} <= end))"
        self._qry = qry

        # Filter to rows for plotting
        plot_gff = self.gff.query(qry).query("feature in @self.gff_features")
        add_gff_fields(plot_gff, ["Parent"])

        # Ensure there are some regions
        assert plot_gff.shape[0] > 0, "No features in this region."

        return plot_gff

    def plot_gff_features(self, 
                          ax, 
                          no_axis=False, 
                          target_gene_id=None
                         ):
        """Plot features of gff in this region"""

        PLUS_STRAND_Y = 2 / 4
        NEG_STRAND_Y = 1 / 4
        
        
        DNA_COLOR = "lightgrey"
        ORF_COLOR = "dimgrey"
        GENE_COLOR = "silver"
        TARGET_COLOR = "teal"

        # Plot features
        for _, row in self.plot_gff.iterrows():

            # Define color and size from feature
            if row["feature"] == "protein_coding_gene":
                lw = 3
                color = GENE_COLOR
            elif row["feature"] == "CDS":
                lw = 8
                color = ORF_COLOR
                if target_gene_id is not None and str(row["Parent"]).startswith(
                    target_gene_id
                ):
                    color = TARGET_COLOR

            # Define y position from strand
            if row["strand"] == "+":
                ypos = PLUS_STRAND_Y
            elif row["strand"] == "-":
                ypos = NEG_STRAND_Y

            # Plot the feature
            ax.plot([row["start"], row["end"]], [ypos, ypos], lw=lw, color=color)

            # Could add annotation text

        # Indicate strands themselves
        ax.plot(
            [self.start, self.end],
            [PLUS_STRAND_Y, PLUS_STRAND_Y],
            lw=1,
            color=DNA_COLOR,
            zorder=-10,
        )
        ax.annotate(
            xy=(self.end, PLUS_STRAND_Y),
            xycoords="data",
            ha="left",
            va="center",
            fontsize=8,
            text=" ($+$)",
        )

        # Indicate strands themselves
        ax.plot(
            [self.start, self.end],
            [NEG_STRAND_Y, NEG_STRAND_Y],
            lw=1,
            color=DNA_COLOR,
            zorder=-10,
        )
        ax.annotate(
            xy=(self.end, NEG_STRAND_Y),
            xycoords="data",
            ha="left",
            va="center",
            fontsize=8,
            text=" ($-$)",
        )
        
        if self.amplicon_region is not None:
            chrom, interval = self.amplicon_region.split(":")
            start, end = interval.split("-")
            start = int(start)
            end = int(end)
            
            ax.plot([start, end],
                    [3/4, 3/4],
                    color='red',
                    lw=4
                   )
            

        # Clean axis ticks
        if no_axis:
            ax.axis("off")
        else:
            for s in ["top", "right", "bottom", "left"]:
                ax.spines[s].set_visible(False)
            ax.get_yaxis().set_ticks([])

        ax.set_xlim((self.start, self.end))
        ax.set_ylim((0, 1))
        ax.label_outer()

        return None


class SequencePlotter:

    # Colors
    AT_PER_COL = "steelblue"
    HP_LENGTH_COL = "firebrick"
    SEQ_COMP_PALETTE = "Greens" 

    def __init__(self, seq):
        """
        Plot various summary statistics of a nucleotide sequennce

        """
        self.seq = seq
        self._calc_summaries()

    def _calc_summaries(self):
        """
        Calculate summary statistics of the sequence

        """

        self.hp_runs = get_homopolymer_runs(self.seq)
        self.per_gc = calc_sliding_percentGC(self.seq, window=20)
        self.seq_array = get_array_encoding(self.seq)

    def plot_sequence_array(self, ax, start, end):
        """
        Plot array giving nucleotide composition

        """
        # Plot
        ax.imshow(
            self.seq_array, cmap=self.SEQ_COMP_PALETTE, aspect="auto", extent=(start, end, 3.5, -0.5)
        )

        # Ticks
        ax.set_yticks([0, 1, 2, 3])
        ax.set_yticklabels(["A", "T", "C", "G"])

        # Labels
        ax.set_xlabel("Position (bp)")
        # ax.label_outer()

        return None

    def plot_sequence_complexity(self, ax, start, end):
        """
        Plot homopolymer run length and AT (%) on a single
        axis

        """

        # Define x values
        xs = np.arange(start, end)

        # Homopolymers
        ax.plot(xs, self.hp_runs, lw=1, color=self.HP_LENGTH_COL, label="Homopolymer Length (bp)")
        ax.set_ylabel("Homopolymer \nLength (bp)", color=self.HP_LENGTH_COL)

        # Limits
        ax.set_ylim((1, ax.get_ylim()[1]))
        ax.set_xlim(start, end)

        # Ticks
        ax.yaxis.set_major_locator(plt.MultipleLocator(10))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(5))

        # GC
        ax.patch.set_visible(False)
        axm = ax.twinx()
        axm.set_zorder(ax.get_zorder() - 1)
        axm.fill_between(
            x=xs,
            y1=0,
            y2=100 * (1 - self.per_gc),
            alpha=0.5,
            color=self.AT_PER_COL,
            label="% AT",
        )
        axm.set_ylim(ax.get_ylim()[0], 100)
        axm.set_ylabel("AT (%)\n[20bp sliding average]", color=self.AT_PER_COL)

        # Ticks
        axm.yaxis.set_major_locator(plt.MultipleLocator(20))
        axm.yaxis.set_minor_locator(plt.MultipleLocator(10))

        # Clean axis
        axm.xaxis.set_visible(False)
        plt.setp(ax.get_xticklabels(), visible=False)


class PrimerPlotter:
    def __init__(self, primer_df):

        # Set primer data frame
        self.primer_df = primer_df

        # Group into pairs, get colors
        self._get_primer_pair_groups()
        self._set_colors()

    def _get_primer_pair_groups(self):
        self.grps = self.primer_df.groupby("pair_name")
        self.n_grps = len(self.grps)
        self.group_names = list(self.grps.groups.keys())

    def _set_colors(self, cmap="Reds"):
        self.col_dt = dict(zip(self.group_names, sns.color_palette(cmap, self.n_grps)))

    def plot(self, ax, start, end, no_axis=False):
        """
        Plot, scaled

        """

        AMP_LW = 3
        PRIMER_LW = 3
        PRIMER_SZ = 40

        for ix, (pair_name, pair_df) in enumerate(self.grps):

            # Extract information
            F_row = pair_df.query("direction == 'F'")
            R_row = pair_df.query("direction == 'R'")

            F_start = F_row["start"]
            F_end = F_start + F_row["length"]

            R_start = R_row["start"]
            R_end = R_start - R_row["length"]

            # Forward
            ax.plot(
                [F_start, F_end], [ix, ix], lw=PRIMER_LW, color=self.col_dt[pair_name]
            )
            ax.scatter(
                x=F_end, y=ix, marker=9, s=PRIMER_SZ, color=self.col_dt[pair_name]
            )

            # Reverse
            ax.plot(
                [R_start, R_end], [ix, ix], lw=PRIMER_LW, color=self.col_dt[pair_name]
            )
            ax.scatter(
                x=R_end, y=ix, marker=8, s=PRIMER_SZ, color=self.col_dt[pair_name]
            )

            # Amplicon
            ax.plot(
                [F_start, R_start], [ix, ix], lw=AMP_LW, color=self.col_dt[pair_name]
            )

            # Annotate
            ax.annotate(
                xy=(R_start, ix),
                ha="left",
                va="center",
                fontsize=6,
                text=f"   {pair_name} {F_row['product_bp'].values[0]}bp",
            )
            # color=col_dt[pair_name])

            # Clean ticks
            ax.set_xlim(start, end)
            if no_axis:
                ax.axis("off")
            else:
                for s in ["top", "right", "bottom", "left"]:
                    ax.spines[s].set_visible(False)
                ax.get_yaxis().set_ticks([])
            ax.label_outer()


class CoverageLengthPlotter:

    # Colors
    READ_LENGTH_PALETTE = "Spectral_r"

    def __init__(self, bam_path):
        """
        Plot a coverage distribution from a `bam_path`

        - Include sanity checks

        """
        self.bam_path = bam_path
        self.bam_df = self._get_bam_dataframe()

    def _get_bam_dataframe(self):
        """
        Create a dataframe from the BAM file

        """
        with pysam.AlignmentFile(self.bam_path, "r") as bam:
            bam_df = pd.DataFrame(
                [
                    AlignmentPlotData.from_pysam_aligned_segment(alignment)
                    for alignment in bam
                ]
            )
        return bam_df

    def _get_positional_dataframe(self):
        """
        Create a dataframe where each row is a position

        """
        dfs = []
        for _, row in self.bam_df.iterrows():
            try:
                position = np.arange(row["ref_start"], row["ref_end"] + 1)
                length = np.repeat(row["query_alignment_length"], position.shape[0])
            except ValueError:
                print("This row FAILED!")
                print(row)
                continue
            dfs.append(pd.DataFrame({"position": position, "length": length}))
        positional_df = pd.concat(dfs)
        return positional_df

    def calc_readlength_summary(self, max_bp=5000, intv_bp=1000):
        """
        Create a dataframe summarising read lengths overlapping
        each position

        """

        positional_df = self._get_positional_dataframe()

        # Bin read lengths
        self.bins = np.arange(0, max_bp + intv_bp, intv_bp)
        positional_df["length_bin"] = pd.cut(
            positional_df["length"], bins=self.bins, right=True, include_lowest=True
        )

        # Summarise read length counts per position
        self.summary_df = (
            positional_df.groupby(["length_bin", "position"]).size().reset_index()
        )
        self.summary_df.rename({0: "count"}, axis=1, inplace=True)
        bin_cats = self.summary_df["length_bin"].dtype.categories[::-1]
        self.summary_df["length_bin"] = pd.Categorical(
            values=self.summary_df["length_bin"],
            categories=bin_cats,
        )
        
        # Define colors from binns
        cols = sns.color_palette(self.READ_LENGTH_PALETTE, len(bin_cats))
        self.bin_col = dict(zip(bin_cats, cols))


    def plot(self, ax, start, end, include_legend=False):
        """
        Plot the coverage by read length histrogram

        """

        # Iterate over length bins, create stacked histogram
        s, e = self.summary_df["position"].min(), self.summary_df["position"].max()
        last_count = np.zeros(e - s + 1)
        for ix, (name, bin_df) in enumerate(self.summary_df.groupby("length_bin")):

            # Extract
            positions = np.array(bin_df["position"])
            counts = np.array(bin_df["count"])
            # assert positions.shape[0] == end - start + 1

            # Plot
            ax.fill_between(
                x=positions,
                y1=last_count,
                y2=last_count + counts,
                zorder=+ix,
                clip_on=False,
                color=self.bin_col[name],
                label=f"{np.abs(name.left):.0f}",
            )

            # Update last count
            last_count += counts

        # Limits
        ax.set_xlim(start, end)
        ax.set_ylim((0, ax.get_ylim()[1]))
        
        # Ticks
        ax.set_axisbelow(True)
        # ax.yaxis.set_major_locator(plt.MultipleLocator(50))
        # ax.yaxis.set_minor_locator(plt.MultipleLocator(10))
        # ax.grid(axis='y', which='major', ls='dotted')

        # Labels
        ax.set_xlabel("Position (bp)")
        ax.set_ylabel("Coverage")

        # Legend
        if include_legend:
            ax.legend(bbox_to_anchor=(1, 1), loc="upper left")
        ax.label_outer()


# --------------------------------------------------------------------------------
# Combined plotting
#
# --------------------------------------------------------------------------------


class CombinedPlotter:
    """Runs, but super ugly"""

    def __init__(self, sequence_plotter, gff_plotter, coverage_plotter, primer_plotter):
        self.seq_plotter = sequence_plotter
        self.gff_plotter = gff_plotter
        self.coverage_plotter = coverage_plotter
        self.primer_plotter = primer_plotter

    def plot(self, start, end, title=None, output_path=None):
        """
        Create a combined plot of sequence composition, complexity,
        gene locations, and candidate primers

        """

        # Constants
        SCALING = 0.08
        # ROWS_PER_PRIMER = 2

        # Axes order and sizing
        AxesFrame = namedtuple("AxesFrame", ["name", "rows", "plot_func"])
        axes_order = [
            AxesFrame(
                "coverage",
                24,
                partial(self.coverage_plotter.plot, start=start, end=end, include_legend=True),
            ),
            # AxesFrame(
            #     "primer", 2, partial(self.primer_plotter.plot, start=start, end=end, no_axis=True)
            # ),
            AxesFrame("genes", 8, partial(self.gff_plotter.plot_gff_features, no_axis=True)),
            AxesFrame(
                "complexity",
                16,
                partial(
                    self.seq_plotter.plot_sequence_complexity, start=start, end=end
                ),
            ),
            AxesFrame(
                "seq",
                8,
                partial(self.seq_plotter.plot_sequence_array, start=start, end=end),
            ),
        ]

        # Figure size
        total_rows = sum([a.rows for a in axes_order]) + len(axes_order)
        height = total_rows * SCALING
        width = 10

        # Create figure
        fig = plt.figure(figsize=(width, height))
        fig.subplots_adjust(hspace=0.2)

        # Create grid
        gs = GridSpec(nrows=total_rows, ncols=1)

        # Iterate over axes and plot
        l = 0
        for i, axis_item in enumerate(axes_order):

            # Prepare axis
            ax = plt.subplot(gs[l : (l + axis_item.rows)])

            # Plot
            axis_item.plot_func(ax)

            # Optionally add title
            if i == 0 and title is not None:
                ax.set_title(title, loc="left")

            # Define boundary of next plot
            l = l + axis_item.rows + 1

        # Optionally write
        if output_path is not None:
            fig.savefig(output_path, bbox_inches="tight", pad_inches=0.5, dpi=300)


# --------------------------------------------------------------------------------
# Main
#
# --------------------------------------------------------------------------------


@click.command(short_help="Plot read coverage over a target gene.")
@click.option(
    "-b",
    "--bam_path",
    default=None,
    required=True,
    help="BAM to use for coverage plot."
    "Should consist of mapped reads ONLY; i.e. from target-extraction."
)
@click.option(
    "-i",
    "--target_id",
    default="PF3D7_0709000",
    show_default=True,
    help="Gene ID of target region.",
)  # or do I need target name from multiplex?
@click.option(
    "-n",
    "--target_name",
    default="CRT1",
    show_default=True,
    help="Gene ID as defined in multiplex.",
)
@click.option(
    "-a",
    "--amplicon_bed",
    default="resources/multiply/multiplexes/multiplex.03.greedy.bed",
    help="BED file defining amplicon positions for multiplex.",
)
@click.option(
    "-p",
    "--primer_csv",
    default="resources/multiply/multiplexes/multiplex.03.greedy.csv",
    help="CSV file defining primer binding sites for multiplex.",
)
@click.option(
    "-r",
    "--reference_name",
    default="Pf3D7",
    show_default=True,
    type=click.Choice(reference_collection),
    help="Reference genome.",
)
@click.option(
    "-w",
    "--window_size",
    default=5000,
    show_default=True,
    help="Size of window to plot, in basepairs.",
)
def main(
    bam_path,
    target_id,
    target_name,
    amplicon_bed,
    primer_csv,
    reference_name,
    window_size,
):
    """
    Create a plot of coverage from a BAM file over a given `target_id`,
    with an amplicon defined by a `multiplex_path`

    """
    # INPUT PARAMETERS
    print("Input parameters")
    print(f"  BAM: {bam_path}")
    print(f"  Target ID: {target_id}")
    print(f"  Target Name: {target_name}")
    print(f"  Amplicon BED: {amplicon_bed}")
    print(f"  Primer CSV: {primer_csv}")
    print(f"  Reference: {reference_name}")
    print(f"  Plotting window size (bp): {window_size}")
    print("Done.\n")

    # Get reference
    reference = reference_collection[reference_name]

    # LOAD GFF
    print("Loading GFF...")
    gff = load_gff(reference.gff_path)
    gff.rename({"seqid": "seqname"}, axis=1, inplace=True)
    print(f"  Number of records: {gff.shape[0]}")
    print("Done.\n")

    # COMPUTE AMPLICON REGION
    print("Preparing plotting region")

    # Get amplicon information from multiplex...
    multiplex_df = pd.read_csv(
        amplicon_bed, names=["seqname", "start", "end", "target_name"], sep="\t"
    )
    assert (
        target_name in multiplex_df["target_name"].tolist()
    ), f"Target {target_name} not in multiplex {multiplex_df['target_name'].tolist()}."
    target_info = multiplex_df.query("target_name == @target_name").squeeze()
    amp_length = target_info["end"] - target_info["start"]
    window_excess = window_size - amp_length
    assert (
        window_size > amp_length
    ), "Window size {window_size}bp, must be larger than amplicon size {amp_length}bp."
    window_start = target_info["start"] - int(window_excess / 2)
    window_end = target_info["end"] + int(window_excess / 2)
    chrom = target_info["seqname"]
    region = f"{chrom}:{window_start}-{window_end-1}"

    # Load sequence data from FASTA...
    seq = load_haplotype_from_fasta(reference.fasta_path, region=region)
    assert window_size == len(
        seq
    ), f"Sequence found is length {len(seq)}bp, but window size was {window_size}."
    print(f"  Amplicon length: {amp_length}")
    print(f"  Region: {region}")
    print(f"  Sequence: {seq[:10]}...")
    print("Done.\n")

    # Get primer information...
    primer_df = pd.read_csv(primer_csv).query("target == @target_id")
    primer_df.insert(3, "pair_name", [s[:-2] for s in primer_df["primer_name"]])
    primer_df.rename({"position": "start", "primer_bp": "length"}, inplace=True, axis=1)

    # INSTANTIATE INDIVIDUAL PLOTTERS
    # Gene format file
    target_region = f"{target_info['seqname']}:{target_info['start']}-{target_info['end']}"
    gff_plotter = GffPlotter(
        gff=gff, 
        chrom=chrom, 
        start=window_start, 
        end=window_end, 
        amplicon_region=target_region)

    # Coverage
    coverage_plotter = CoverageLengthPlotter(bam_path=bam_path)
    coverage_plotter.calc_readlength_summary()

    # Sequence plotter
    sequence_plotter = SequencePlotter(seq)

    # Primer plotter
    primer_plotter = PrimerPlotter(primer_df=primer_df)

    # CREATE COMBINED PLOT
    combined_plotter = CombinedPlotter(
        sequence_plotter=sequence_plotter,
        gff_plotter=gff_plotter,
        coverage_plotter=coverage_plotter,
        primer_plotter=primer_plotter,
    )
    combined_plotter.plot(
        start=window_start,
        end=window_end,
        title=f"{target_id} | {target_name}",
        output_path=f"scripts/figures/diagram.{target_id}.{target_name}.pdf",
    )


if __name__ == "__main__":
    main()
