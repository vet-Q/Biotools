import os
import click
import subprocess
import pandas as pd
from collections import namedtuple
from nomadic.lib.generic import produce_dir, print_header, print_footer


# ================================================================
# Run MAFFT
#
# ================================================================


def run_mafft(input_fasta, output_msa, method="mafft"):
    """Run MAFFT on an input fasta file"""

    assert method in ["mafft", "linsi", "ginsi"]
    cmd = f"{method} {input_fasta} > {output_msa}"
    subprocess.run(cmd, shell=True, check=True)

    return None


def load_msa_as_dict(msa_path):
    """
    Load MSA from MAFFT as a dictionary

    """

    dt = {}
    with open(msa_path, "r") as msa:
        line = msa.readline()
        while line:
            if line.startswith(">"):
                header = line.strip()
                dt[header] = ""
            else:
                dt[header] += line.strip()
            line = msa.readline()

    # Sanity check
    assert all([len(v) for k, v in dt.items()])

    return dt


# ================================================================
# Call SNPs from a MSA
#
# ================================================================


def create_snp_table(ref_msa, samp_msa, add_chrom=None):
    """
    Given two sequences from a multiple sequence alignment,
    produce a table of SNPs

    """

    # Ensure length equal
    assert len(ref_msa) == len(samp_msa), "MSAs should be same length."

    # Define record
    Record = namedtuple("Record", ["POS", "ID", "REF", "ALT"])

    # Iterate
    idx = 0
    records = []
    for r, s in zip(ref_msa, samp_msa):
        if r == "-":
            continue
        if s == "-":
            idx += 1
            continue
        if r == s:
            idx += 1
            continue

        # Record SNP
        records.append(Record(idx, ".", r, s))
        idx += 1
    snp_table = pd.DataFrame(records)

    # Optionally add a chromosome indicator
    if add_chrom:
        snp_table.insert(0, "CHROM", add_chrom)

    return snp_table


def build_joint_snp_table_from_msa(msa_path):
    """
    Build a combined variant table for all alignments
    within an `msa_path`

    """

    # Load the MSA as a dictionary
    dt = load_msa_as_dict(msa_path)

    # Extract reference strain information
    ref_header = [h for h in dt if h.startswith(">Pf3D7")][0]
    ref_msa = dt.pop(ref_header)
    ref_region = ref_header.split("|")[-1].strip()
    ref_contig, ref_intv = ref_region.split(":")
    ref_start, ref_end = ref_intv.split("-")

    # For each alignment in the MSA, create a variant table
    var_dfs = []
    for header, msa in dt.items():

        # Parse information
        strain = header.split("|")[0][1:].strip()
        region = header.split("|")[-1].strip()

        # Get a variant data frame
        var_df = create_snp_table(ref_msa, msa, add_chrom=ref_contig)

        if len(var_df) > 0:

            var_df.insert(0, "strain", strain)
            var_df.insert(1, "region", region)

            # Adjust to reference coordinates
            print(var_df.columns)
            var_df["POS"] += int(ref_start)

            # Store
            var_dfs.append(var_df)

    if var_dfs:
        return pd.concat(var_dfs)

    else:
        print("No variants discovered.")
        return pd.DataFrame([])


# ================================================================
# Main script
#
# ================================================================


@click.command(short_help="Create an MSA and call SNPs.")
@click.option(
    "-f",
    "--fasta_dir",
    type=click.Path(exists=True),
    default="resources/truthsets/fastas",
    help="Path to directory containing .fasta files.",
)
def msacall(fasta_dir):
    """
    Run a MSA (via MAFFT), and then call SNPs with respect to the
    reference genome

    """
    # Define directories
    t0 = print_header("TRUTHSET: Create an MSA and call SNPs")
    msa_dir = produce_dir(fasta_dir.replace("fasta", "msa"))

    fastas = [f"{fasta_dir}/{fasta}" for fasta in os.listdir(fasta_dir) if fasta.endswith(".fasta")]
    print(f"Found {len(fastas)} .fasta files in {fasta_dir}.")

    # RUN MAFFT
    print("Creating MSAs with MAFFT...")
    output_msas = []
    for fasta in fastas:
        output_msa = fasta.replace(".fasta", ".mafft.aln").replace("fasta", "msa")
        run_mafft(fasta, output_msa)
        output_msas.append(output_msa)
    print("Done.")
    print("")

    # CALL SNPs
    snp_dir = produce_dir(msa_dir.replace("msas", "msa_snps"))
    print("Creating genotype CSVs...")
    for input_msa in output_msas:

        print(f"  {input_msa}")

        # Find SNPs
        target_gt_df = build_joint_snp_table_from_msa(input_msa)

        # Write, if any SNPs were found
        if target_gt_df.shape[0] > 0:
            target_gt_df.to_csv(
                f"{snp_dir}/{os.path.basename(input_msa).replace('.mafft.aln', '.snp.csv')}"
            )
    print("Done.")
    print("")
    print_footer(t0)
