import os
import numpy as np
import pandas as pd
import pysam
import subprocess

from dataclasses import dataclass

from nomadic.lib.process_gffs import load_gff
from nomadic.lib.process_gffs import add_gff_fields

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

plt.rcParams["figure.dpi"] = 100


# ================================================================================
# PARAMETERS
#
# ================================================================================


EXPT_DIR = "../experiments/2021-11-14_strain-validation-flongle-lfb"
FOCUS_BARCODE = 1

MULTIPLEX_PATH = "../resources/multiply/multiplexes/multiplex.03.greedy.csv"
GFF_PATH = "../resources/plasmodb/52/PlasmoDB-52_Pfalciparum3D7.gff"

WINDOW_SIZE_BP = 5000
N_READS = 30

OUTPUT_DIR = "figures/fig1_visualise-reads"
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)


# ================================================================================
# LOADING
#
# ================================================================================


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


def load_bam_info(bam_path):
    with pysam.AlignmentFile(bam_path, "r") as bam:
        bam_df = pd.DataFrame(
            [
                AlignmentPlotData.from_pysam_aligned_segment(alignment)
                for alignment in bam
            ]
        )
    return bam_df


# ================================================================================
# PLOTTING
#
# ================================================================================


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
        add_gff_fields(plot_gff, ["Parent"])

        # Ensure there are some regions
        assert plot_gff.shape[0] > 0, "No features in this region."

        return plot_gff

    def plot_gff_features(self, ax, no_axis=False, target_gene_id=None):
        """Plot features of gff in this region"""

        PLUS_STRAND_Y = 3 / 4
        NEG_STRAND_Y = 1 / 4

        # Plot features
        for _, row in self.plot_gff.iterrows():

            # Define color and size from feature
            if row["feature"] == "protein_coding_gene":
                lw = 3
                color = "darkgrey"
            elif row["feature"] == "CDS":
                lw = 7
                color = "teal"
                if target_gene_id is not None and str(row["Parent"]).startswith(
                    target_gene_id
                ):
                    color = "darkorange"

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
            color="lightgrey",
            zorder=-10,
        )
        ax.annotate(
            xy=(self.end, PLUS_STRAND_Y),
            xycoords="data",
            ha="left",
            va="center",
            fontsize=6,
            text=" ($+$)",
        )

        # Indicate strands themselves
        ax.plot(
            [self.start, self.end],
            [NEG_STRAND_Y, NEG_STRAND_Y],
            lw=1,
            color="lightgrey",
            zorder=-10,
        )
        ax.annotate(
            xy=(self.end, NEG_STRAND_Y),
            xycoords="data",
            ha="left",
            va="center",
            fontsize=6,
            text=" ($-$)",
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


def offset_tick_formatter(position_bp, num_kbp_ix):
    """
    Plot as kbp offset from first tick
    
    """
    if num_kbp_ix == 0:
        return int(position_bp)
    return f"+{int(num_kbp_ix)}kbp"


# ================================================================================
# MAIN
#
# ================================================================================


def main():
    """
    Create visualisations of a random set of reads in genomic
    context for all amplicons from a multiplex

    """

    # Set experiment directories
    nomadic_dir = f"{EXPT_DIR}/nomadic/guppy/hac/single_end/"
    targets_dir = f"{nomadic_dir}/barcodes/barcode{FOCUS_BARCODE:02d}/target-extraction"
    print(f"Focus barcode: {FOCUS_BARCODE}")

    # Load GFF file
    gff = load_gff(GFF_PATH)
    gff.rename({"seqid": "seqname"}, axis=1, inplace=True)

    # Load multiplex dataframe
    multiplex_df = pd.read_csv(MULTIPLEX_PATH)
    target_ids = np.unique(multiplex_df["target"])

    # Iterate over target gene IDs
    for target_id in target_ids:
        # Query for target
        target_df = multiplex_df.query("target == @target_id")

        # Extract properties
        target_name = target_df.query("direction == 'F'")["gene_name"].iloc[0]    
        if target_name == "MDR1part":
            target_name = "MDR1"
        start = int(target_df.query("direction == 'F'")["position"])
        end = int(target_df.query("direction == 'R'")["position"])
        amplicon_size_bp = end - start
        chrom_int = int(target_id.split("_")[1][:2])
        chrom = f"Pf3D7_{chrom_int:02d}_v3"

        print(f"Target: {target_id} | {target_name}")
        print(f"  Start: {start}")
        print(f"  End: {end}")
        print(f"  Amplicon length (bp): {amplicon_size_bp}")
        print(f"  Chromosome: {chrom}\n")

        # Compute pad sizes
        pad_size = (WINDOW_SIZE_BP - amplicon_size_bp) / 2

        # Load BAM file
        bam_path = f"{targets_dir}/reads.target.{target_name}.bam"
        #subprocess.run(f"tabix {bam_path}", shell=True, check=True) already indexed
        bam_df = load_bam_info(bam_path)
        print(f"Loaded {bam_df.shape[0]} reads.")

        # Subsample and sort reads
        read_df = bam_df.sample(N_READS)
        read_df.sort_values(["query_alignment_length"], ascending=False, inplace=True)
        read_df.insert(0, "plot_ix", range(N_READS))
        

        # Instantiate GFF plot class
        gff_plotter = GffPlotter(
            gff=gff, start=start - pad_size, end=end + pad_size, chrom=chrom
        )

        # Plot
        print("Plotting...")
        output_path = f"{OUTPUT_DIR}/reads.plot.{target_name}.pdf"

        # CREATE FIGURE
        fig = plt.figure(figsize=(7, 2))
        fig.subplots_adjust(hspace=0.05)

        gs = GridSpec(nrows=8, ncols=1)
        ax_gff = plt.subplot(gs[-1])
        ax_reads = plt.subplot(gs[:-1], sharex=ax_gff)

        # READ PLOT
        for _, row in read_df.iterrows():
            ax_reads.plot(
                [row["ref_start"], row["ref_end"]],
                [row["plot_ix"], row["plot_ix"]],
                color="darkgrey",
            )
            ax_reads.axis("off")

        # Prepare title
        title = f"{target_name} | {target_id} | {amplicon_size_bp}bp"
        ax_reads.set_title(title, loc="left")

        # GFF PLOT
        gff_plotter.plot_gff_features(ax_gff, no_axis=False, target_gene_id=target_id)
        ax_gff.xaxis.set_major_locator(plt.LinearLocator(int(WINDOW_SIZE_BP/1000) + 1))
        ax_gff.xaxis.set_major_formatter(plt.FuncFormatter(offset_tick_formatter))

        if output_path is not None:
            fig.savefig(output_path, bbox_inches="tight", pad_inches=0.1)

        print("Done.\n")


if __name__ == "__main__":
    main()
