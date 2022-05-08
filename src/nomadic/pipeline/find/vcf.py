import pysam
import pandas as pd
import allel
from collections import namedtuple


def load_vcf_as_df(vcf_path, only_snps=True):
    """Load a VCF file and return as pandas DataFrame"""

    # Define information to extract
    # - I want to add genotyping information here
    # - Probably also (maybe) want allelic depth information
    Record = namedtuple("VCFRecord", ["chrom", "pos", "ref", "alt", "qual"])

    # Extract
    with pysam.VariantFile(vcf_path) as vcf:

        records = []
        for r in vcf:

            if only_snps:

                if "INDEL" in r.info:
                    if r.info["INDEL"]:
                        continue

                if len(r.ref) > 1:
                    continue

            # Store
            records.append(
                Record(r.chrom, int(r.pos), r.ref.upper(), r.alts[0].upper(), r.qual)
            )

    return pd.DataFrame(records)


def load_vcf_using_allel(vcf_path, add_genotypes=True, only_snp=True):
    """
    Load VCF file for a *single* sample using scikit-allel,
    optionally adding genotypes

    TODO:
    - Handle indels or non-bialellic alleles

    """

    # Load callset
    callset = allel.read_vcf(vcf_path)

    if callset is None:
        print(f"No variants found in file at {vcf_path}.")
        columns = ["chrom", "pos", "ref", "alt", "filter", "qual"]
        columns += ["gt"] if add_genotypes else []
        return pd.DataFrame(columns=columns)

    # Check number of samples
    assert callset["samples"].shape[0] == 1, "Function built for 1 sample."

    # Convert to dataframe
    vcf_df = pd.DataFrame(
        {
            "chrom": callset["variants/CHROM"],
            "pos": callset["variants/POS"],
            "ref": callset["variants/REF"],
            "alt": callset["variants/ALT"][:, 0],
            "filter": callset["variants/FILTER_PASS"],
            "qual": callset["variants/QUAL"],
        }
    )

    # Add genotypes
    if add_genotypes:
        gts = [gt.sum() for gt in callset["calldata/GT"][:, 0]]
        vcf_df["gt"] = gts

    # Remove any indels
    if only_snp:
        ref_snp = [len(r) == 1 for r in vcf_df["ref"]]
        alt_snp = [len(a) == 1 for a in vcf_df["alt"]]
        vcf_df = vcf_df[[r and a for r, a in zip(ref_snp, alt_snp)]]

    return vcf_df
