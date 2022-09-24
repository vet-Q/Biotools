import click
import re
import os
import numpy as np
import pandas as pd

from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.process_gffs import load_gff, add_gff_fields
from nomadic.lib.references import (
    PlasmodiumFalciparum3D7,
    PlasmodiumFalciparumDd2,
    PlasmodiumFalciparumGB4,
    PlasmodiumFalciparumHB3,
)
from nomadic.truthset.fasta import load_haplotype_from_fasta


# ================================================================
# Parameters
#
# ================================================================

REFERENCE_COLLECTION = {
    r.name: r
    for r in [
        PlasmodiumFalciparum3D7(),
        PlasmodiumFalciparumDd2(),
        PlasmodiumFalciparumGB4(),
        PlasmodiumFalciparumHB3(),
    ]
}
REFERENCE_NAME = "Pf3D7"
DEFAULT_CSV_PATH = "resources/multiply/multiplexes/multiplex.03.greedy.csv"
DEFAULT_OUTPUT_DIR = "resources/truthsets/stratifications"


# ================================================================
# Functions
#
# ================================================================


def write_gff_to_bed(df, bed_path, bed_columns=["seqid", "start", "end", "ID"]):
    """Convert from GFF to a bed"""
    sep = "\t"
    with open(bed_path, "w") as bed:
        for _, row in df.iterrows():
            bed.write(f"{sep.join([str(v) for v in row[bed_columns]])}\n")


def get_homopolymer_encoding(seq, l_max=10):
    """
    For a given sequence `seq`, produce a homopolymer
    block size encoding
    i.e.
    ATTCCC = 1, 2, 2, 3, 3, 3
    """

    # Store
    L = len(seq)
    h = np.zeros(L, "int8")
    h[:] = -1

    # Initialise
    i = 0
    p = seq[i]
    l = 0

    # Iterate
    for c in seq:
        if c == p:
            l += 1
        else:
            h[i : (i + l)] = l
            i += l
            p = c
            l = 1

    # Terminate
    h[i : (i + l)] = l
    h[h > l_max] = l_max
    assert (h > 0).all(), "Error in homopolymer encoding."
    assert (h <= l_max).all(), "Error in homopolymer encoding."

    return h


# ================================================================
# Main script
#
# ================================================================


@click.command(short_help="Create stratifications for happy.")
@click.option(
    "-c",
    "--csv_path",
    type=click.Path(exists=True),
    default=DEFAULT_CSV_PATH,
    help="Path to multiplex primer CSV.",
)
def stratify(csv_path):
    """
    Create stratifications of an amplicon panel for
    hap.py comparison analysis

    Steps


    """

    # PARSE INPUTS
    script_descrip = "TRUTHSET: Create stratifications for hap.py analysis"
    t0 = print_header(script_descrip)

    # CREATE OUTPUT DIRECTORY
    print("Inputs")
    output_dir = produce_dir(DEFAULT_OUTPUT_DIR)
    print(f"  Output directory: {output_dir}")
    print(f"  Multiplex primer CSV: {csv_path}")
    print(f"  Reference: {REFERENCE_NAME}")
    print("Done.")
    print("")

    # PREPARE STORAGE
    bed_files = {}

    # LOAD MULTIPLEX
    print("Loading multiplex CSV...")
    csv_name = os.path.basename(csv_path)
    multiplex_df = pd.read_csv(csv_path)
    target_ids = multiplex_df["target"].unique()
    print(f"Found {len(target_ids)} targets...")
    print("Done.")
    print("")

    # CREATE MULTIPLEX BED
    print("Creating amplicon BED file...")
    amplicon_bed_path = (
        f"{output_dir}/{csv_name.replace('.csv', '.full_amplicons.bed')}"
    )
    amplicon_regions = {}
    with open(amplicon_bed_path, "w") as bed:
        for target_id, target_df in multiplex_df.groupby("target"):
            # Chromsome
            pattern = "PF3D7_([0-9]){2}[0-9]+"
            chrom_int = int(re.match(pattern, target_id).group(1))
            chrom = f"Pf3D7_{chrom_int:02d}_v3"  # NOT GENERAL!

            # Start and end
            start = target_df.query("direction == 'F'")["position"].values[0] - 1
            end = target_df.query("direction == 'R'")["position"].values[0]
            assert start < end

            # Write
            line = "\t".join([chrom, str(start), str(end), target_id])
            bed.write(f"{line}\n")

            # Store amplicon regions
            amplicon_regions[target_id] = f"{chrom}:{start}-{end}"
    bed_files["amplicon"] = amplicon_bed_path
    print("Done.")
    print("")

    # LOAD GFF
    print("Loading GFF...")
    reference = REFERENCE_COLLECTION[REFERENCE_NAME]
    gff_df = load_gff(gff_path=reference.gff_path)
    gff_df = add_gff_fields(gff_df)
    print("Done.")
    print("")

    # CREATE CDS BEDs
    print("Creating CDS-based BED files...")
    # Filter to CDS entries
    cds_df = gff_df.query("feature == 'CDS'")
    cds_df["Parent"] = [p.split(".")[0] for p in cds_df["Parent"]]  # clean
    target_cds_df = cds_df.query("Parent in @target_ids")
    print(
        f"  Found coding sequences for {len(target_cds_df['Parent'].unique())} targets."
    )
    print(f"  Totalling {target_cds_df.shape[0]} exons.")
    print("Writing...")

    # Write
    detailed_cds_bed_path = amplicon_bed_path.replace(
        "full_amplicons.bed", "cds_detailed.bed"
    )
    write_gff_to_bed(
        df=target_cds_df,
        bed_path=detailed_cds_bed_path,
        bed_columns=["seqid", "start", "end", "ID"],
    )
    bed_files["exons"] = detailed_cds_bed_path
    cds_bed_path = amplicon_bed_path.replace("full_amplicons.bed", "cds.bed")
    write_gff_to_bed(
        df=target_cds_df,
        bed_path=cds_bed_path,
        bed_columns=["seqid", "start", "end", "Parent"],
    )
    bed_files["cds"] = cds_bed_path
    print(f"  By-target breakdown: {cds_bed_path}")
    print(f"  By-exon breakdown: {detailed_cds_bed_path}")
    print("Done.")
    print("")

    # GET SEQUENCES
    # - Iterate over regions
    # - Compute HP encoding
    # - Convert to BED
    # print("Loading sequences for homopolymer stratifications...")
    # for target_id, region in amplicon_regions.items():
    #     seq = load_haplotype_from_fasta(
    #         fasta_path=reference.fasta_path,
    #         region=region
    #     )
    #     print(f"  Target: {target_id}")
    #     print(f"  Sequence: {seq[:10]}...")

    #     # Calculate homopolymers
    #     print("  Computing homopolymer encoding...")
    #     hp = get_homopolymer_encoding(seq)

    # WRITE STRATIFICATIONS FILE
    stratifications_path = f"{output_dir}/stratification.tsv"
    with open(stratifications_path, "w") as file:
        for bn, bf in bed_files.items():
            file.write(f"{bn}\t{os.path.basename(bf)}\n")
    print(f"Stratificaitons file: {stratifications_path}")

    print_footer(t0)
