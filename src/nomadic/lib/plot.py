from nomadic.lib.statistics import (
    get_homopolymer_runs,
    calc_sliding_percentGC,
    get_array_encoding,
)
import numpy as np
import matplotlib.pyplot as plt


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
        #ax.label_outer()

        return None

    def plot_sequence_complexity(self, ax, start, end):
        """
        Plot homopolymer run length and AT (%) on a single
        axis

        """

        # Define x values
        xs = np.arange(start, end)

        # ----------------------------------------
        # A+T (%) PLOT
        ax.fill_between(
            x=xs,
            y1=0,
            y2=100 * (1 - self.per_gc),
            alpha=0.5,
            color=self.AT_PER_COL,
            label="% AT",
        )
        # Limits
        ax.set_xlim(start, end)
        ax.set_ylim(0, 100)
        # Labels
        lab = "AT (%)"
        #lab += "\n[20bp sliding average]"
        ax.set_ylabel("AT (%)", color=self.AT_PER_COL)
        # Ticks
        ax.yaxis.set_major_locator(plt.MultipleLocator(20))
        ax.yaxis.set_minor_locator(plt.MultipleLocator(10))

        # ----------------------------------------
        # Homopolymer Plot
        # Twin
        axm = ax.twinx()
        # Plot
        axm.plot(xs, self.hp_runs, lw=1, color=self.HP_LENGTH_COL, label="Homopolymer Length (bp)")
        # Limits
        axm.set_xlim(start, end)
        axm.set_ylim(1, 45)
        # Labels
        axm.set_ylabel("Homopolymer \nLength (bp)", color=self.HP_LENGTH_COL)
        # Ticks
        axm.yaxis.set_major_locator(plt.MultipleLocator(10))
        axm.yaxis.set_minor_locator(plt.MultipleLocator(5))

        # Clean axis
        axm.xaxis.set_visible(False)
        plt.setp(ax.get_xticklabels(), visible=False)

     

        # # Homopolymers
        # ax.plot(xs, self.hp_runs, lw=1, color=self.HP_LENGTH_COL, label="Homopolymer Length (bp)")
        # ax.set_ylabel("Homopolymer \nLength (bp)", color=self.HP_LENGTH_COL)

        # # Limits
        # ax.set_ylim((1, ax.get_ylim()[1]))
        # ax.set_xlim(start, end)

        # # GC
        # ax.patch.set_visible(False)
        # axm = ax.twinx()
        # axm.set_zorder(ax.get_zorder() - 1)
        # axm.fill_between(
        #     x=xs,
        #     y1=0,
        #     y2=100 * (1 - self.per_gc),
        #     alpha=0.5,
        #     color=self.AT_PER_COL,
        #     label="% AT",
        # )
        # axm.set_ylim(ax.get_ylim()[0], 100)
        # axm.set_ylabel("AT (%)\n[20bp sliding average]", color=self.AT_PER_COL)

        # # Clean axis
        # axm.xaxis.set_visible(False)
        # plt.setp(ax.get_xticklabels(), visible=False)


class GffPlotter:

    gff_features = ["protein_coding_gene", "CDS"]

    def __init__(self, gff, chrom, start, end):
        """
        Plot information from a GFF over a defined
        region

        """
        self.gff = gff
        self.chrom = chrom
        self.start = start
        self.end = end

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

        # Ensure there are some regions
        assert plot_gff.shape[0] > 0, "No features in this region."

        return plot_gff

    def plot_gff_features(self, ax):
        """Plot features of gff in this region"""
        
        PLUS_STRAND_Y = 3/4
        NEG_STRAND_Y = 1/4

        # Plot features
        for _, row in self.plot_gff.iterrows():

            # Define color and size from feature
            if row["feature"] == "protein_coding_gene":
                lw = 3
                color = "darkgrey"
            elif row["feature"] == "CDS":
                lw = 8
                color = "teal"

            # Define y position from strand
            if row["strand"] == "+":
                ypos = PLUS_STRAND_Y
            elif row["strand"] == "-":
                ypos = NEG_STRAND_Y

            # Plot the feature
            ax.plot([row["start"], row["end"]], [ypos, ypos], lw=lw, color=color)

            # Could add annotation text
            
        # Indicate strands themselves
        ax.plot([self.start, self.end], 
                [PLUS_STRAND_Y, PLUS_STRAND_Y],
                lw=1, color='lightgrey', zorder=-10)
        ax.annotate(xy=(self.end, PLUS_STRAND_Y), 
                    xycoords="data",
                    ha="left", va="center",
                    fontsize=8,
                    text=" ($+$) strand")
        
        # Indicate strands themselves
        ax.plot([self.start, self.end], 
                [NEG_STRAND_Y, NEG_STRAND_Y],
                lw=1, color='lightgrey', zorder=-10)
        ax.annotate(xy=(self.end, NEG_STRAND_Y), 
                    xycoords="data",
                    ha="left", va="center",
                    fontsize=8,
                    text=" ($-$) strand")

        # Clean axis ticks
        ax.axis("off")
        ax.set_xlim((self.start, self.end))
        ax.set_ylim((0, 1))
        ax.label_outer()

        return None