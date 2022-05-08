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
    """
    Implementation for calling variants with `longshot`

    """

    def set_arguments(self, fasta_path):
        """
        Set any required arguments as instance variables

        """
        self.fasta_path = fasta_path

    def call_variants(self):
        """
        Call variants

        """
        cmd = "longshot -F"
        cmd += f" --bam {self.bam_path}"
        cmd += f" --ref {self.fasta_path}"
        cmd += f" --out {self.vcf_path}"

        subprocess.run(cmd, shell=True, check=True)


class BcfTools(VariantCaller):
    """
    Implementation for calling variants with `bcftools`

    """

    def set_arguments(self, fasta_path):
        """
        Set any required arguments as instance variables

        """
        self.fasta_path = fasta_path

    def call_variants(self):
        """
        Call variants

        """
        cmd = "bcftools mpileup -Ou"
        cmd += f" -f {self.fasta_path}"
        cmd += f" {self.bam_path}"
        cmd += " | bcftools call -cv"  # consensus calling
        cmd += f" -Ou -o {self.vcf_path}"

        subprocess.run(cmd, shell=True, check=True)


class GatkHaplotypeCaller(VariantCaller):
    """
    Implementation for calling variants with `gatk HaplotypeCaller`;
    for now, just running with default settings

    """

    with_module_load = True  # at present, only running from cluster
    module = "GATK/4.2.5.0-GCCcore-11.2.0-Java-11"

    def set_arguments(self, fasta_path):
        """
        Set any required arguments as instance variables

        """
        self.fasta_path = fasta_path

    def call_variants(self):
        """
        Call variants

        """
        cmd = ""
        if self.with_module_load:
            cmd = f"module load {self.module} && "

        cmd += "gatk HaplotypeCaller"
        cmd += f" -R {self.fasta_path}"
        cmd += f" -I {self.bam_path}"
        cmd += f" -O {self.vcf_path}"

        subprocess.run(cmd, shell=True, check=True)


# ================================================================
# Define collection of available callers
#
# ================================================================

# Note, they are already initialised
caller_collection = {
    "longshot": LongShot(),
    "bcftools": BcfTools(),
    "gatk": GatkHaplotypeCaller(),
}
