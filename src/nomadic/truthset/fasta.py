import click
import json
import subprocess
import pysam
import pandas as pd
from nomadic.lib.generic import produce_dir, print_header, print_footer
from nomadic.lib.references import (
    PlasmodiumFalciparum3D7,
    PlasmodiumFalciparumDd2,
    PlasmodiumFalciparumGB4,
    PlasmodiumFalciparumHB3
)

# ================================================================
# Parameters
#
# ================================================================


PAD_BP = 8_000
CSV_PATH = "resources/truthsets/table.target_genes.csv"


# ================================================================
# Loading functions
#
# ================================================================


def load_haplotype_from_fasta(fasta_path, **kwargs):
    """
    Load a particular region of a fasta file

    """

    if kwargs is None:
        fasta_seq = next(pysam.FastxFile(fasta_path)).sequence
    else:
        fasta_seq = pysam.FastaFile(fasta_path).fetch(**kwargs)

    return fasta_seq


@click.command(short_help="Generate FASTA files for targets.")
@click.option(
    "-g",
    "--gff_path",
    type=str,
    default=CSV_PATH,
    help="Path to .gff file of *only* targets."
)
@click.option(
    "-p",
    "--pad_bp",
    type=int,
    default=PAD_BP,
    help="Amount of additional sequencing to pad around target for in FASTA."
)
def fasta(gff_path, pad_bp):
    """
    Create FASTA files for every target in a `gff_path`, across all strains
    
    """

    # Output directory
    t0 = print_header("TRUTHSET: Produce FASTA files for targets")
    output_dir = produce_dir("resources/truthsets", "fastas")
    print(f"Output directory: {output_dir}")
    print(f"Target .gff path: {gff_path}")
    print(f"Pad size (bp): {pad_bp}")
    print("")

    # Define reference genomes
    references = [
        PlasmodiumFalciparum3D7(),
        PlasmodiumFalciparumDd2(),
        PlasmodiumFalciparumGB4(),
        PlasmodiumFalciparumHB3()
    ]
    references_dt = {r.name: r for r in references}

    # Load combined target data frame
    combined_df = pd.read_csv(gff_path)
    assert "standard_id" in combined_df.columns, "Must have a 'standard_id' column."

    # Create a FASTA for each target
    for target_id, target_df in combined_df.groupby("standard_id"):

        target_names = ", ".join(list(target_df["name"].unique()))

        print("."*80)
        print(f"Target: {target_id} | {target_names}")
        print("."*80)

        # Open fasta
        target_fasta_path = f"{output_dir}/{target_id}.fasta"
        print(f"Output: {target_fasta_path}")
        with open(target_fasta_path, "w") as fasta:

            for _, row in target_df.iterrows():

                # Extract reference
                reference = references_dt[row["reference"]]
                print(f"  Writing FASTA for {reference.name}...")

                # Define header
                annot = f"{row['reference']} | {row['ID']} | {row['name']} | "
                region = f"{row['seqid']}:{row['start']-1-pad_bp}-{row['end']+pad_bp}"
                header = f">{annot}{region}"

                # Extract sequence
                seq = load_haplotype_from_fasta(fasta_path=reference.fasta_path, region=region)

                # To stdout
                print(f"  Header: {header}")
                print(f"  Seq: {seq[:10]} ... {seq[-10:]} = {len(seq)}bp")
                print("")

                # Write to FASTA
                lines = f"{header}\n{seq}\n"
                fasta.write(lines)
        print("Done.")
        print("")
    print_footer(t0)