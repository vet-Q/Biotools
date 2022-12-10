import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


from nomadic.lib.references import (
    PlasmodiumFalciparum3D7,
    PlasmodiumFalciparumDd2,
    PlasmodiumFalciparumGB4,
    PlasmodiumFalciparumHB3
)
from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.generic import produce_dir, print_header, print_footer
from nomadic.pipeline.coi.trim.targets import TARGET_COLLECTION


# Load PAF
def load_paf(paf_path: str) -> pd.DataFrame:
    """Load Pairwise Alignment File"""
    
    # Columns
    PAF_COLUMNS = [
        "query_name",
        "query_length",
        "query_start",
        "query_end",
        "strand_match",
        "target_name",
        "target_length",
        "target_start",
        "target_end",
        "n_matches_alignment",
        "n_bases_alignment",
        "mapq"
    ]
    
    # Dataframe
    paf_df = pd.read_csv(paf_path, header=None, sep="\t")
    paf_df = paf_df.iloc[:, :12]
    paf_df.columns = PAF_COLUMNS
    paf_df["identity"] = paf_df["n_matches_alignment"] / paf_df["n_bases_alignment"]
    
    return paf_df


# From to Plasmodb
EXPECTED_SIZES = {
    "MSP2": {
        "Pf3D7": 272 * 3,
        "PfDd2": 296 * 3,
        "PfGB4": 292 * 3,
        "PfHB3": 256 * 3
    }
}

def plot_target_histogram(read_df, 
                          target_gene, 
                          references,
                          palette="Set1",
                          plot_expected=True,
                          output_path=None):
    """
    Plot read length histograms for a target gene,
    annotate by best mapping reference
    
    """
    
    assert "highest_identity_ref" in read_df.columns
    assert "length" in read_df.columns
    
    # CREATE COLOR PALETTE
    ref_cols = dict(
        zip([r.name for r in references], 
            sns.color_palette(palette, len(references)))
    )
    
    # GROUP BY REFERENCE
    grps = read_df.groupby("highest_identity_ref")
    
    # PLOT
    fig, ax = plt.subplots(1, 1, figsize=(6, 3.5))

    ax.hist(
        [gdf["length"] for _, gdf in grps],
        bins=50,
        stacked=True,
        color=[ref_cols[r] for r, _ in grps],
        ec="black",
        lw=0.5,
        label=[r for r, _ in grps]
    )

    # Labels
    ax.set_xlabel(f"{target_gene} ORF length (bp)")
    ax.set_ylabel("Count")

    # Ticks and grid
    ax.set_axisbelow(True)
    ax.xaxis.set_major_locator(plt.MultipleLocator(50))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(10))
    ax.grid(ls='dotted', alpha=0.5)

    plot_expected=True
    if plot_expected and target_gene in EXPECTED_SIZES:
        for ref, l in EXPECTED_SIZES[target_gene].items():
            ax.axvline(l, color=ref_cols[ref], lw=2, zorder=-1)
            ax.annotate(xy=(l, ax.get_ylim()[1]),
                        ha="center", va="bottom",
                        rotation=90,
                        fontsize=8,
                        color=ref_cols[ref],
                        text=f" {l}bp"
                       )

    # Legend
    ax.legend(title="Highest Identity", frameon=False)

    # Save
    if output_path is not None:
        fig.savefig(output_path, dpi=300, pad_inches=0.5, bbox_inches="tight")


def main(expt_dir, config, barcode, target_gene):
    """
    Filter and trim all reads in a BAM file to overlap a `target_gene` 
    and span `start` and `end` positions; then convert to FASTQ
    
    """

    # PARSE INPUTS
    script_descrip = "NOMADIC: Trim target BAM file for COI analysis"
    t0 = print_header(script_descrip)
    script_dir = "coi"
    params = build_parameter_dict(expt_dir, config, barcode)

    target = TARGET_COLLECTION[target_gene]
    print("User inputs:")
    print(f"  Target: {target.name}")
    print(f"  Chrom: {target.chrom}")
    print(f"  Start: {target.start}")
    print(f"  End: {target.end}")
    print("Done.\n")

    # Focus on a single barcode, if specified
    if "focus_barcode" in params:
        params["barcodes"] = [params["focus_barcode"]]
        
    # References
    references = [
        PlasmodiumFalciparum3D7(),
        PlasmodiumFalciparumDd2(),
        PlasmodiumFalciparumGB4(),
        PlasmodiumFalciparumHB3()
    ]

    # ITERATE
    print("Iterating over barcodes...")
    for barcode in params["barcodes"]:
        print("." * 80)
        print(f"Barcode: {barcode}")
        print("." * 80)

        # DIRECTORIES
        # Inputs
        barcode_dir = f"{params['barcodes_dir']}/{barcode}"
        coi_dir = f"{barcode_dir}/coi"

        # Output
        plot_dir = produce_dir(coi_dir, "plots")

        # Clipped read information
        fastq_csv_path = f"{coi_dir}/fastq_clipped/reads.target.{target_gene}.clipped.csv"
        read_df = pd.read_csv(fastq_csv_path)

        # Load Panmap dataframes
        panmap_dfs = []
        for r in references:
            input_paf = f"{coi_dir}/panmap/{barcode}.{r.name}.{target_gene}.sorted.paf"
            paf_df = load_paf(input_paf)
            paf_df.insert(0, "barcode", barcode)
            paf_df.insert(1, "reference", r.name)
            panmap_dfs.append(paf_df)
        panmap_df = pd.concat(panmap_dfs)

        # Merge in highest identity across panel mapping
        panmap_wide_df = pd.pivot(
            index="query_name",
            columns="reference",
            values="identity",
            data=panmap_df
        )
        highest_identity_ref = panmap_wide_df.idxmax(axis=1)
        highest_identity_ref.name = "highest_identity_ref"

        read_df = pd.merge(left=read_df, 
                        right=highest_identity_ref,
                        left_on="read_id",
                        right_index=True
                        )

        # Plot read length histogram
        plot_target_histogram(
            read_df, 
            target_gene, 
            references,
            palette="Set1",
            output_path=f"{plot_dir}/plot.read_lengths.{target_gene}.pdf")

    print_footer(t0)




