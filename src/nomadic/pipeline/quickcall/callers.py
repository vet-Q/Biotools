import os
import shutil
import subprocess
from abc import ABC, abstractmethod
from nomadic.lib.process_vcfs import bcftools_reheader, bcftools_index


# ================================================================
# Define abstract base class for different variant calling methods
#
# ================================================================


class VariantCaller(ABC):
    def __init__(self, fasta_path: str) -> None:
        self.fasta_path = fasta_path

    @abstractmethod
    def _run(self, bam_path: str, vcf_path: str) -> None:
        pass

    def run(self, bam_path: str, vcf_path: str, sample_name: str = None):
        """
        Run core variant calling method

        Also:
        - Optionally adding sample name
        - Indexing the output VCF

        """

        # Store
        self.vcf_path = vcf_path

        # Core method
        self._run(bam_path, vcf_path)

        # Optionally name
        if sample_name is not None:
            bcftools_reheader(self.vcf_path, self.vcf_path, [sample_name])

        # Index
        bcftools_index(self.vcf_path)


# ================================================================
# Concrete implementations
#
# ================================================================


class BcfTools(VariantCaller):
    # SETTINGS
    ANNOTATE = "FORMAT/DP,FORMAT/AD"
    MAX_DEPTH = 10_000
    MIN_DEPTH = 40
    MIN_QUAL = 20

    def _run(self, bam_path: str, vcf_path: str) -> None:
        """
        Run variant calling with bcftools

        TODO:
        Note that here I am following recommendations of
        `bcftools` in calling with ONT reads, that is, using
        `-X ont` for `bcftools mpileup`, which sets:

        `-B`  : disable per-base alignment quality
        `-Q5` : Skip bases with base quality < 5
        `--max-BQ 30` : Set the maximum base quality to 30
            - ONT sets homopolymers to 90 for some reason
        `-I`  : Skip indel calling
            - I only report SNPs in dashboard
            - But could be difficult for some variants (K76T)

        Then for `bcftools call`, I use `-P 0.01`,

        `-P` : Prior on mutation rate

        """

        cmd_pileup = "bcftools mpileup -Ou"
        cmd_pileup += f" -X ont"  # set to ONT mode.
        cmd_pileup += f" --annotate {self.ANNOTATE}"
        cmd_pileup += f" --max-depth {self.MAX_DEPTH}"
        cmd_pileup += f" -f {self.fasta_path}"
        cmd_pileup += f" {bam_path}"

        cmd_call = "bcftools call -mv -P 0.01 - "

        cmd_filter = "bcftools view"
        cmd_filter += " --min-alleles 2"
        cmd_filter += " --max-alleles 2"
        cmd_filter += " --types='snps'"
        cmd_filter += f" -e 'FORMAT/DP<{self.MIN_DEPTH}||QUAL<{self.MIN_QUAL}' -"

        cmd_tags = f"bcftools +fill-tags -Oz -o {vcf_path} - -- -t FORMAT/VAF"

        cmd = f"{cmd_pileup} | {cmd_call} | {cmd_filter} | {cmd_tags}"

        subprocess.run(cmd, check=True, shell=True)


# ================================================================
# Define collection of available callers
#
# ================================================================


# Note, they are already initialised
caller_collection = {
    "bcftools": BcfTools,
}
