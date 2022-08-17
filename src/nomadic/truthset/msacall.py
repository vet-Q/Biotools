import os
import click
import subprocess
import pandas as pd
from functools import reduce
from collections import namedtuple
from nomadic.lib.generic import produce_dir, print_header, print_footer
from nomadic.lib.process_vcfs import (
    bcftools_view,
    bcftools_sort,
    bcftools_index,
    bcftools_concat,
    bcftools_query_samples,
)
from nomadic.lib.references import PlasmodiumFalciparum3D7


REFERENCE = PlasmodiumFalciparum3D7()

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


# def create_snp_table(ref_msa, samp_msa, add_chrom=None):
#     """
#     Given two sequences from a multiple sequence alignment,
#     produce a table of SNPs

#     """

#     # Ensure length equal
#     assert len(ref_msa) == len(samp_msa), "MSAs should be same length."

#     # Define record
#     Record = namedtuple("Record", ["POS", "ID", "REF", "ALT"])

#     # Iterate
#     idx = 0
#     records = []
#     for r, s in zip(ref_msa, samp_msa):
#         if r == "-":
#             continue
#         if s == "-":
#             idx += 1
#             continue
#         if r == s:
#             idx += 1
#             continue

#         # Record SNP
#         records.append(Record(idx, ".", r, s))
#         idx += 1
#     snp_table = pd.DataFrame(records)

#     # Optionally add a chromosome indicator
#     if add_chrom:
#         snp_table.insert(0, "CHROM", add_chrom)

#     return snp_table


# def build_joint_snp_table_from_msa(msa_path, reference_name="Pf3D7"):
#     """
#     Build a combined variant table for all alignments
#     within an `msa_path`

#     """

#     # Load the MSA as a dictionary
#     dt = load_msa_as_dict(msa_path)

#     # Extract reference strain information
#     ref_header = [h for h in dt if h.startswith(reference_name)][0]
#     ref_msa = dt.pop(ref_header)
#     ref_region = ref_header.split("|")[-1].strip()
#     ref_contig, ref_intv = ref_region.split(":")
#     ref_start, ref_end = ref_intv.split("-")

#     # For each alignment in the MSA, create a variant table
#     var_dfs = []
#     for header, msa in dt.items():

#         # Parse information
#         strain = header.split("|")[0][1:].strip()
#         region = header.split("|")[-1].strip()

#         # Get a variant data frame
#         var_df = create_snp_table(ref_msa, msa, add_chrom=ref_contig)

#         if len(var_df) > 0:

#             var_df.insert(0, "strain", strain)
#             var_df.insert(1, "region", region)

#             # Adjust to reference coordinates
#             print(var_df.columns)
#             var_df["POS"] += int(ref_start)

#             # Store
#             var_dfs.append(var_df)

#     if var_dfs:
#         return pd.concat(var_dfs)

#     else:
#         print("No variants discovered.")
#         return pd.DataFrame([])


# ================================================================
# Facilitate MSA to VCF conversion
#
# ================================================================


def create_snp_df(ref_seq, samp_seq):
    """
    Given two sequences from a multiple sequence alignment,
    produce a table of SNPs

    """

    # Ensure length equal
    assert len(ref_seq) == len(samp_seq), "MSAs should be same length."

    # Define record
    columns = ["POS", "ID", "REF", "ALT"]
    Record = namedtuple("Record", columns)

    # Iterate
    idx = 0
    records = []
    for r, s in zip(ref_seq, samp_seq):
        if r == "-":
            continue
        if s == "-":
            idx += 1
            continue
        if r == s:
            idx += 1
            continue

        # Record SNP
        records.append(Record(int(idx), ".", r, s))
        idx += 1
    snp_table = pd.DataFrame(records, columns=columns)

    return snp_table


class MSAtoVCF:

    vcf_columns = [
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
    ]
    vcf_sep = "\t"

    def __init__(self, msa_path):
        """
        Convert a multiple sequence alignment into a VCF file

        """
        self.msa_path = msa_path
        self.msa_dt = load_msa_as_dict(msa_path)

    def set_reference(self, reference_name="Pf3D7"):
        """
        Indicate which of the sequences in the MSA is the reference genome

        """

        # Assign reference name
        self.reference_name = reference_name

        # Get reference sequence
        self.reference_seq = [
            seq for header, seq in self.msa_dt.items() if self.reference_name in header
        ]
        if not len(self.reference_seq) == 1:
            raise ValueError(
                f"Exactly one sequence in MSA must include {reference_name} in header, found {len(self.reference_seq)}"
            )
        self.reference_seq = self.reference_seq[0]

        # Extract reference metadata (used to build VCF)
        self.reference_header = [
            header for header in self.msa_dt if self.reference_name in header
        ][0]
        self.reference_region = self.reference_header.split("|")[-1].strip()
        self.reference_contig, self.reference_intv = self.reference_region.split(":")
        self.reference_start, self.reference_end = self.reference_intv.split("-")

    def _create_snp_dfs(self):
        """
        For every pairwise sample, find all the snps

        """

        self.snp_dfs = []
        for header, seq in self.msa_dt.items():
            snp_df = create_snp_df(ref_seq=self.reference_seq, samp_seq=seq)
            sample_name = header.split("|")[0][1:].strip()  # assumption about header
            snp_df[sample_name] = "1/1"  # SNPs that are found are homozygous ALT
            self.snp_dfs.append(snp_df)

    def _merge_snp_dfs(self, fill_missing="0/0"):
        """
        Merge all of the SNP data frames into a single table,
        filling missing entires with homozygous reference


        """
        # Merge all data frames in the list
        self.merged_snp_df = reduce(
            lambda df1, df2: pd.merge(
                left=df1, right=df2, on=["POS", "ID", "REF", "ALT"], how="outer"
            ),
            self.snp_dfs,
        ).fillna(fill_missing)

    def _format_as_vcf(self):
        """
        Format the merged SNP dataframe into a minimal VCF

        """

        # Convert bases to upper case
        self.merged_snp_df["REF"] = [c.upper() for c in self.merged_snp_df["REF"]]
        self.merged_snp_df["ALT"] = [c.upper() for c in self.merged_snp_df["ALT"]]

        # Split out variant information
        variant_df = self.merged_snp_df[
            [c for c in self.merged_snp_df.columns if c in self.vcf_columns]
        ]

        # Add required VCF columns
        variant_df["CHROM"] = self.reference_contig
        variant_df["POS"] = variant_df["POS"] + int(self.reference_start)
        variant_df["QUAL"] = "."
        variant_df["FILTER"] = "PASS"
        variant_df["INFO"] = "."
        variant_df["FORMAT"] = "GT"
        variant_df = variant_df[self.vcf_columns]

        # Split out sample information
        sample_df = self.merged_snp_df[
            [c for c in self.merged_snp_df.columns if c not in self.vcf_columns]
        ]

        # Create VCF dataframe
        self.vcf_df = pd.concat([variant_df, sample_df], axis=1)

    @staticmethod
    def _create_vcf_header(reference):
        """
        Create a valid header for a VCF file,
        including `contig` lines for the `reference`

        """

        # Essential format specification
        header = "##fileformat=VCFv4.2\n"

        with open(reference.fasta_path, "r") as fasta:
            for contig in fasta:

                # Skip if not header
                if not contig.startswith(">"):
                    continue

                # Parse header
                ID, _, _, length, _ = contig.split("|")
                ID = ID[1:].strip()
                length = int(length.split("=")[1].strip())

                # Add to string
                header += f"##contig=<ID={ID},length={length}>\n"

        # Genotype fields
        header += "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>\n"

        return header

    def create_vcf(self, output_path=None):
        """
        Create a VCF file from a multiple sequence alignment (MSA),
        after a reference has been set

        NB: Need a valid header, including genotype fields and contigs,
        in order for downstream bcftools commands (e.g. bcftools sort),
        to work.

        """
        assert self.reference_name is not None, "Please `.set_reference()` first."

        self._create_snp_dfs()
        self._merge_snp_dfs()
        self._format_as_vcf()

        # Optionally write the VCF file
        if output_path is not None:
            with open(output_path, "w") as vcf:

                # Write the header
                vcf.write(self._create_vcf_header(REFERENCE))

                # Write the information fields
                vcf.write(f"#{self.vcf_sep.join(self.vcf_df.columns)}\n")
                for _, row in self.vcf_df.iterrows():
                    vcf.write(f"{self.vcf_sep.join([str(v) for v in row.values])}\n")


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
    output_dir = produce_dir("resources", "truthsets", "mafft")
    msa_dir = produce_dir(output_dir, "msas")

    fastas = [fasta for fasta in os.listdir(fasta_dir) if fasta.endswith(".fasta")]
    print(f"Found {len(fastas)} .fasta files in {fasta_dir}.")

    # RUN MAFFT
    print("Creating MSAs with MAFFT...")
    output_msas = []
    for fasta in fastas:
        input_fasta = f"{fasta_dir}/{fasta}"
        output_msa = f"{msa_dir}/{fasta.replace('.fasta', '.mafft.aln')}"
        run_mafft(input_fasta, output_msa)
        output_msas.append(output_msa)
    print("Done.")
    print("")

    # CALL SNPs
    vcf_dir = produce_dir(msa_dir.replace("msas", "vcfs"))
    vcf_target_dir = produce_dir(vcf_dir, "by_target")
    print("Creating genotype VCFs...")
    vcfs_to_concat = []
    for input_msa in output_msas:

        print(f"  {input_msa}")

        # Write, if any SNPs were found
        vcf_builder = MSAtoVCF(input_msa)
        vcf_builder.set_reference()

        # Writing VCF
        vcf_path = f"{vcf_target_dir}/{os.path.basename(input_msa).replace('.mafft.aln', '.snp.vcf')}"
        print(f"  VCF written to: {vcf_path}")
        vcf_builder.create_vcf(output_path=vcf_path)

        # Sort and index
        sorted_vcf = f"{vcf_path}.sorted.gz"
        bcftools_sort(input_vcf=vcf_path, O="z", output_vcf=sorted_vcf)
        bcftools_index(sorted_vcf)

        # Store
        vcfs_to_concat.append(sorted_vcf)

    # Concatenate and index across genes
    target_genes_vcf = f"{vcf_dir}/target_genes.concat.vcf.gz"
    bcftools_concat(input_vcfs=vcfs_to_concat, output_vcf=target_genes_vcf)
    bcftools_index(target_genes_vcf)

    # Split to individual sample VCFs
    samples = bcftools_query_samples(input_vcf=target_genes_vcf)
    for sample in samples:
        vcf_sample_path = target_genes_vcf.replace(".concat", f".{sample}.concat")
        bcftools_view(
            input_vcf=target_genes_vcf, output_vcf=vcf_sample_path, s=sample, O="z", dry_run=False
        )
        bcftools_index(vcf_sample_path)

    print(f"Multi-sample concatenated VCF written to: {target_genes_vcf}")

    print("Done.")
    print("")
    print_footer(t0)
