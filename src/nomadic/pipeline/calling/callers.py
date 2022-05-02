import subprocess
from abc import ABC, abstractmethod


# ================================================================
# Define abstract base class for different variant calling methods
#
# ================================================================


class VariantCaller(ABC):
    """
    Base class for a selection of variant callers

    """

    def __init__(self):
        self.bam_path = None
        self.vcf_path = None

    def set_files(self, bam_path, vcf_path):
        self.bam_path = bam_path
        self.vcf_path = vcf_path

    @abstractmethod
    def set_arguments(self):
        pass

    @abstractmethod
    def call_variants(self):
        pass


# ================================================================
# Concrete implementations
#
# ================================================================


class LongShot(VariantCaller):
    def set_arguments(self, fasta_path):
        """
        Set command line arguments necessary to run
        longshot

        """

        self.fasta_path = fasta_path

        return None

    def call_variants(self):
        """
        Call variants using longshot

        """

        cmd = "longshot -F"
        cmd += f" --bam {self.bam_path}"
        cmd += f" --ref {self.fasta_path}"
        cmd += f" --out {self.vcf_path}"

        subprocess.run(cmd, shell=True, check=True)

        return None


class BcfTools(VariantCaller):
    def set_arguments(self, fasta_path):
        """Set bcftools call arguments"""

        self.fasta_path = fasta_path

        return None

    def call_variants(self):
        """Call variants using bcftools call"""

        cmd = "bcftools mpileup -Ou"
        cmd += f" -f {self.fasta_path}"
        cmd += f" {self.bam_path}"
        cmd += " | bcftools call -cv"  # consensus calling
        cmd += f" -Ou -o {self.vcf_path}"

        subprocess.run(cmd, shell=True, check=True)

        return None


# ================================================================
# Define collection of available callers
#
# ================================================================

# Note, they are already initialised
caller_collection = {"longshot": LongShot(), "bcftools": BcfTools()}
