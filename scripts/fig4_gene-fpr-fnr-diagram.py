import pandas as pd

from nomadic.lib.generic import produce_dir
from nomadic.lib.references import reference_collection
from nomadic.lib.plot import SequencePlotter, GffPlotter
from nomadic.lib.process_gffs import load_gff, add_gff_fields
from nomadic.truthset.fasta import load_haplotype_from_fasta

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


# ================================================================================
# Amplicon limits
#
# ================================================================================


def limit_to_target(target_name, amplicon_df, vcf_df):
    """
    Limit `vcf_df` to rows within a target amplicon

    """

    # Extract relevant target info
    target_info = amplicon_df.query("target_name == @target_name").squeeze()
    chrom = target_info["seqname"]
    start = target_info["start"]
    end = target_info["end"]

    return vcf_df.query("CHROM == @chrom").query("(@start <= POS) and (POS <= @end)")


# ================================================================================
# Plot TPR and FPR
#
# ================================================================================


def plot_tpr_fpr(df, ax):
    """
    Plot calling performance of all True Positive sites,
    via their TPR;
    and those False Positive sites with a FPR > 0

    """

    ax.scatter(
        x="POS",
        y="TPR",
        s="n_reads",
        c="None",
        ec="green",
        lw=0.5,
        data=df.query("TPR > 0"),
        clip_on=False,
    )

    ax.scatter(
        x="POS",
        y="FPR",
        s="n_reads",
        c="None",
        ec="red",
        lw=0.5,
        data=df.query("FPR > 0"),
        clip_on=False,
    )


# ================================================================================
# Main
#
# ================================================================================


def main(combined_vcf_path, output_dir, method, amplicon_bed, window_size=None):
    """
    Create a positional plot of TPR / FPR for all genes in an amplicon
    panel defined by `amplicon_bed`

    """

    # Load combined VCF
    combined_vcf_df = pd.read_csv(combined_vcf_path)

    # COMPUTE FPR / TPR
    n_barcodes = combined_vcf_df["barcode"].unique().shape[0]
    n_replicates = combined_vcf_df["rep"].unique().shape[0]
    m = n_barcodes * n_replicates
    print(n_barcodes, n_replicates, m)

    summary_df = (
        combined_vcf_df.groupby(["n_reads", "CHROM", "POS", "REF", "ALT"])
        .agg(
            N=pd.NamedAgg("BS", len),
            TPC=pd.NamedAgg("TP", sum),
            FPC=pd.NamedAgg("FP", sum),
            FNC=pd.NamedAgg("FN", sum),
        )
        .reset_index()
    )

    summary_df["TPR"] = summary_df["TPC"] / summary_df["N"]
    summary_df["FNR"] = summary_df["FNC"] / summary_df["N"]
    summary_df["FPR"] = summary_df["FPC"] / m

    # LOAD REFERENCE
    reference = reference_collection["Pf3D7"]

    # LOAD GFF
    print("Loading GFF...")
    gff = load_gff(reference.gff_path)
    gff.rename({"seqid": "seqname"}, axis=1, inplace=True)
    print(f"  Number of records: {gff.shape[0]}")
    print("Done.\n")

    # LOAD ALL AMPLICONS
    amplicon_df = pd.read_csv(
        amplicon_bed, names=["seqname", "start", "end", "target_name"], sep="\t"
    )
    amplicon_df.index = amplicon_df["target_name"]

    # ITERATE OVER AMPLICONS
    for target_name in amplicon_df["target_name"]:

        print(f"Plotting for target: {target_name}...")

        # EXTRACT AMPLICON INFORMATION
        seqname, start, end, _ = amplicon_df.loc[target_name]
        amp_length = end - start

        # Optionally compute window size
        if window_size is not None:
            window_excess = window_size - amp_length
            assert window_excess > 0
        else:
            window_excess = 0

        window_start = start - int(window_excess / 2)
        window_end = end + int(window_excess / 2)

        # Define region
        region = f"{seqname}:{window_start}-{window_end-1}"

        # LOAD SEQUENCE
        seq = load_haplotype_from_fasta(reference.fasta_path, region=region)
        print(f"  Amplicon length: {amp_length}")
        print(f"  Region: {region}")
        print(f"  Sequence: {seq[:10]}...")

        # LIMIT TO TARGET MUTATIONS
        target_df = limit_to_target(target_name, amplicon_df, summary_df)
        print(f"  Found {target_df.shape[0]} calls within {target_name}.")
        print(
            f"  Representing {len(target_df.groupby(['CHROM', 'POS']))} unique sites."
        )

        # PLOT
        # --------------------------------------------------------------------------------
        # GFF plotter
        gff_plotter = GffPlotter(
            gff=gff, chrom=seqname, start=window_start, end=window_end
        )

        # Sequence plotter
        seq_plotter = SequencePlotter(seq)

        fig = plt.figure(figsize=(10, 4))

        gs = GridSpec(nrows=9, ncols=1)
        ax_gff = plt.subplot(gs[0])
        ax_tpr = plt.subplot(gs[1:5], sharex=ax_gff)
        ax_hp = plt.subplot(gs[5:7], sharex=ax_gff)
        ax_seq = plt.subplot(gs[7:9], sharex=ax_gff)

        # Plot TPR
        plot_tpr_fpr(df=target_df, ax=ax_tpr)
        ax_tpr.tick_params(labelbottom=False)
        ax_tpr.set_ylim((0, 1))
        ax_tpr.yaxis.set_major_locator(plt.MultipleLocator(0.2))
        ax_tpr.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
        ax_tpr.set_axisbelow(True)
        ax_tpr.grid(ls="dotted")

        # GFF PLOT
        gff_plotter.plot_gff_features(ax=ax_gff)
        ax_gff.set_title(target_name, loc="left")
        ax_gff.tick_params(axis="x", labelbottom=False)

        # HP PLOT
        seq_plotter.plot_sequence_complexity(
            ax=ax_hp, start=window_start, end=window_end
        )

        # SEQ PLOT
        seq_plotter.plot_sequence_array(ax=ax_seq, start=window_start, end=window_end)

        fig.savefig(
            f"{output_dir}/plot.fpr_tpr.{target_name}.{method}.pdf",
            dpi=300,
            bbox_inches="tight",
            pad_inches=0.5,
        )
        # --------------------------------------------------------------------------------
        print("Done.\n")


if __name__ == "__main__":

    # Settings
    combined_vcf_path = "summary_tables/2021-11-14_strain-validation-flongle-lfb/sup/fig4_aggregated-vcfs-clair3sing.csv"
    output_dir = produce_dir(
        "scripts", "figures", "fig4", "2021-11-14_strain-validation-flongle-lfb", "sup"
    )

    # Run
    main(
        combined_vcf_path,
        output_dir,
        method="clair3sing",
        amplicon_bed="resources/multiply/multiplexes/multiplex.03.greedy.bed",
        window_size=4000,
    )
