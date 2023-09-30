import os
import shutil
import subprocess
from abc import ABC, abstractmethod
from nomadic.lib.process_vcfs import bcftools_reheader


# ================================================================
# Define abstract base class for different variant calling methods
#
# ================================================================


class VariantCaller(ABC):
    def __init__(self):
        """
        Abstract base class defining interface for a collection
        of variant calling algorithms

        Facilitates BAM --> VCF

        """
        self.bam_path = None
        self.vcf_path = None

    def set_files(self, bam_path, vcf_path):
        self.bam_path = bam_path
        self.vcf_path = vcf_path

    @abstractmethod
    def set_arguments(self):
        pass

    @abstractmethod
    def call_variants(self, sample_name=None):
        """
        Call variants, optionally naming sample in
        output VCF file

        """
        pass


# ================================================================
# Concrete implementations
#
# ================================================================


class LongShot(VariantCaller):
    """
    Encapsulate `longshot`

    """

    def set_arguments(self, fasta_path):
        """
        Set any required arguments as instance variables

        """
        self.fasta_path = fasta_path

    def call_variants(self, sample_name=None):
        """
        Call variants with `longshot`

        """
        cmd = "longshot -F"
        cmd += f" --bam {self.bam_path}"
        cmd += f" --ref {self.fasta_path}"
        cmd += f" --out {self.vcf_path}"

        if sample_name is not None:
            cmd += f" --sample_id {sample_name}"

        subprocess.run(cmd, shell=True, check=True)


class BcfTools(VariantCaller):
    """
    Encapsulate `bcftools call`; using consensus method (`-cv`)

    """

    MAX_DEPTH = 10_000  # we don't really want to limit this
    ANNOTATE = "FORMAT/DP,FORMAT/AD"

    def set_arguments(self, fasta_path):
        """
        Set any required arguments as instance variables

        """
        self.fasta_path = fasta_path

    def call_variants(self, sample_name=None):
        """
        Call variants with `bcftools call`

        """
        cmd = "bcftools mpileup -Ou"
        cmd += f" --annotate {self.ANNOTATE}"
        cmd += f" --max-depth {self.MAX_DEPTH}"
        cmd += f" -f {self.fasta_path}"
        cmd += f" {self.bam_path}"
        cmd += " | bcftools call -cv"  # consensus calling
        cmd += f" -Ou -o {self.vcf_path}"

        subprocess.run(cmd, shell=True, check=True)

        if sample_name is not None:
            bcftools_reheader(self.vcf_path, self.vcf_path, [sample_name])


class GatkHaplotypeCaller(VariantCaller):
    """
    Encapsulate `gatk HaplotypeCaller`;
    for now, just running with default settings

    """

    with_module_load = True  # at present, only running from cluster
    module = "GATK/4.2.5.0-GCCcore-11.2.0-Java-11"

    def set_arguments(self, fasta_path):
        """
        Set any required arguments as instance variables

        """
        self.fasta_path = fasta_path

    def call_variants(self, sample_name=None):
        """
        Call variants with `gatk HaplotypeCaller`

        """
        cmd = ""
        if self.with_module_load:
            cmd = f"module load {self.module} && "

        cmd += "gatk HaplotypeCaller"
        cmd += f" -R {self.fasta_path}"
        cmd += f" -I {self.bam_path}"
        cmd += f" -O {self.vcf_path}"

        subprocess.run(cmd, shell=True, check=True)


class Clair3Singularity(VariantCaller):
    """
    Encapsulate Clair3; run via a singularity container

    NB:
    - We don't have the exact Guppy model
    - This is specific for Slurm

    """

    SIF_PATH = "/u/jash/containers/clair3_latest.sif"
    BIND_DIRS = True
    THREADS = "4"

    # Guppy models
    # MODEL = "/u/jash/projects/rerio/clair3_models/r1041_e82_400bps_sup_g615" # SUP
    # MODEL = "/u/jash/projects/rerio/clair3_models/r1041_e82_400bps_hac_g632"  # HAC
    
    # Dorado models
    MODEL = "/u/jash/projects/rerio/clair3_models/r1041_e82_400bps_sup_v420"
    # MODEL = "/u/jash/projects/rerio/clair3_models/r1041_e82_400bps_fast_v420"

    def set_arguments(self, fasta_path):
        """
        Set arguments for Clair3

        Note that all paths must be absolute; and any directories that
        do not exist should be created before binding.

        """

        # reference genome
        self.fasta_path = os.path.abspath(fasta_path)

        # paths must be absolute
        self.bam_path = os.path.abspath(self.bam_path)
        self.vcf_path = os.path.abspath(self.vcf_path)

        # clair3 outputs to a directory, not a file
        self.vcf_dir = self.vcf_path.replace(".vcf", "")

        # Create if doesn't exist, helps with mounting (-B)
        if not os.path.exists(self.vcf_dir):
            os.makedirs(self.vcf_dir)

        # Collect directories
        # NB: We need to include the model as well, if we are not using
        # one within the container
        self.dirs = [
            os.path.dirname(self.bam_path),
            os.path.dirname(self.fasta_path),
            self.vcf_dir,
            self.MODEL,
        ]

    def call_variants(self, sample_name=None):
        """
        Call variants with Clair3

        """

        # Build command
        cmd = f"singularity exec"
        if self.BIND_DIRS:
            cmd += f" {' '.join([f'-B {d}' for d in self.dirs])}"
        cmd += f" {self.SIF_PATH} /opt/bin/run_clair3.sh"
        cmd += f" --bam_fn={self.bam_path}"
        cmd += f" --ref_fn={self.fasta_path}"
        cmd += f" --threads={self.THREADS}"
        cmd += " --platform='ont'"
        # cmd += f" --model_path=/opt/models/{self.MODEL}"
        cmd += f" --model_path={self.MODEL}"
        cmd += f" --output {self.vcf_dir}"
        cmd += " --include_all_ctgs"
        cmd += " --enable_phasing"

        if sample_name is not None:
            cmd += f" --sample_name={sample_name}"

        # Run
        subprocess.run(cmd, check=True, shell=True)

        # The next two steps are for consistency with other callers
        # Decompress output
        cmd = f"bcftools view {self.vcf_dir}/phased_merge_output.vcf.gz -Ou -o {self.vcf_dir}/phased_merge_output.vcf"
        subprocess.run(cmd, check=True, shell=True)

        # Move VCF file
        shutil.copyfile(f"{self.vcf_dir}/phased_merge_output.vcf", self.vcf_path)
        shutil.copyfile(
            f"{self.vcf_dir}/phased_merge_output.vcf.gz.tbi", f"{self.vcf_path}.tbi"
        )


class DeepVariantSingularity(VariantCaller):
    """
    Call variants with PEPPER-Margin-DeepVariant via singularity

    """

    SIF_PATH = "tools/pepper_deepvariant_r0.8.sif"
    BIND_DIRS = True
    THREADS = "4"
    MODEL = "ont_r9_guppy5_sup"  # All that is available

    def set_arguments(self, fasta_path):
        """
        Set arguments for PEPPER-Margin-DeepVariant

        Note that all paths must be absolute; and any directories that
        do not exist should be created before binding.

        """

        # reference genome
        self.fasta_path = os.path.abspath(fasta_path)

        # paths must be absolute
        self.bam_path = os.path.abspath(self.bam_path)
        self.vcf_path = os.path.abspath(self.vcf_path)

        # clair3 outputs to a directory, not a file
        self.vcf_dir = self.vcf_path.replace(".vcf", "")
        self.vcf_prefix = os.path.basename(self.vcf_dir)

        # Create if doesn't exist, helps with mounting (-B)
        if not os.path.exists(self.vcf_dir):
            os.makedirs(self.vcf_dir)

        # Collect directories
        self.dirs = [
            os.path.dirname(self.bam_path),
            os.path.dirname(self.fasta_path),
            self.vcf_dir,
        ]

    def call_variants(self, sample_name=None):
        """
        Call variants with PEPPER-Margin-DeepVariant

        """

        # Build command
        cmd = f"singularity exec"
        if self.BIND_DIRS:
            cmd += f" {' '.join([f'-B {d}' for d in self.dirs])}"
        cmd += f" {self.SIF_PATH} run_pepper_margin_deepvariant call_variant"
        cmd += f" -b {self.bam_path}"
        cmd += f" -f {self.fasta_path}"
        cmd += f" -t {self.THREADS}"
        cmd += f" -o {self.vcf_dir}"
        cmd += f" -p {self.vcf_prefix}"
        cmd += f" --{self.MODEL}"

        if sample_name is not None:
            cmd += f" -s {sample_name}"

        # Run
        subprocess.run(cmd, check=True, shell=True)

        # The next two step is for consistency with other callers
        # Decompress output
        cmd = f"bcftools view {self.vcf_dir}/{self.vcf_prefix}.vcf.gz -Ou -o {self.vcf_path}"
        subprocess.run(cmd, check=True, shell=True)

        # Move VCF index file
        shutil.copyfile(
            f"{self.vcf_dir}/{self.vcf_prefix}.vcf.gz.tbi", f"{self.vcf_path}.tbi"
        )


# ================================================================
# Define collection of available callers
#
# ================================================================

# Note, they are already initialised
caller_collection = {
    "longshot": LongShot(),
    "bcftools": BcfTools(),
    "gatk": GatkHaplotypeCaller(),
    "clair3sing": Clair3Singularity(),
    "deepvarsing": DeepVariantSingularity(),
}
