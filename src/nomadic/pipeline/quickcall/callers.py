import os
import shutil
import subprocess
from abc import ABC, abstractmethod
from nomadic.lib.generic import produce_dir
from nomadic.lib.process_vcfs import bcftools_reheader, bcftools_index


# ================================================================
# Define abstract base class for different variant calling methods
#
# ================================================================


class VariantCaller(ABC):
    def __init__(self, fasta_path: str) -> None:
        self.fasta_path = fasta_path
        self.vcf_path = None

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

    def filter(self, 
               output_vcf: str, 
               bed_path: str,
               min_depth: int = 50, 
               #min_qual: int = 20,
               to_biallelic: bool = False
               ) -> None:
        """
        Filters the output VCF to only regions contained within
        `bed_path`; also optionally excludes sites below a threshold

        """

        if self.vcf_path is None:
            raise ValueError("Must run variant calling before filtering.")

        self.MIN_DEPTH = min_depth
        #self.MIN_QUAL = min_qual

        cmd_view = "bcftools view"
        cmd_view += f" -R {bed_path}"
        if to_biallelic:
            cmd_view += " --types='snps'"
            cmd_view += " --min-alleles 2"
            cmd_view += " --max-alleles 2"
        cmd_view += f" {self.vcf_path}"

        cmd_filter = "bcftools filter"
        cmd_filter += " -S ."
        cmd_filter += f" -e 'FORMAT/DP<{self.MIN_DEPTH}'" # ||QUAL<{self.MIN_QUAL}'"
        cmd_filter += f" -Oz -o {output_vcf} -"

        cmd = f"{cmd_view} | {cmd_filter}"

        subprocess.run(cmd, check=True, shell=True)

        # Index
        bcftools_index(output_vcf)


# ================================================================
# Concrete implementations
#
# ================================================================


class BcfTools(VariantCaller):
    # SETTINGS
    ANNOTATE_MPILEUP = "FORMAT/DP,FORMAT/AD"
    ANNOTATE_CALL = "FORMAT/GQ"
    MAX_DEPTH = 10_000

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
        cmd_pileup += f" --annotate {self.ANNOTATE_MPILEUP}"
        cmd_pileup += f" --max-depth {self.MAX_DEPTH}"
        cmd_pileup += f" -f {self.fasta_path}"
        cmd_pileup += f" {bam_path}"

        # NB: We are returning *all* variants (not using -v)
        cmd_call = f"bcftools call -a 'FORMAT/GQ' -m -P 0.01 -Oz -o {vcf_path} -"
        # #cmd_call += f" --annotate {self.ANNOTATE_CALL}"
        # cmd_call += " -Oz -o {vcf_path} -"

        cmd = f"{cmd_pileup} | {cmd_call}"

        subprocess.run(cmd, check=True, shell=True)


class Clair3Singularity(VariantCaller):
    """
    Run Clair3 using Singularity

    Note that Singularity only runs natively on Linux

    """

    SIF_PATH = "/u/jash/containers/clair3_latest.sif"

    # In theory, the model should match the version of the basecalling
    # software (guppy/dorado) that was used
    MODEL = "/u/jash/projects/rerio/clair3_models/r1041_e82_400bps_sup_v420"

    def _mount_dirs(self, bam_path: str, vcf_path: str) -> None:
        """
        Prepare all directories and path names for mounting

        """
        self.bam_path = os.path.abspath(bam_path)
        self.vcf_path = os.path.abspath(vcf_path)
        self.fasta_path = os.path.abspath(self.fasta_path)

        # Clair3 outputs to a directory, not a file
        # - create if doesn't exist, helps with mounting (-B)
        self.vcf_dir = produce_dir(self.vcf_path.replace(".vcf.gz", ""))

        # Collect directories
        # NB: We need to include the model as well, if we are not using
        # one within the container
        self.dirs = [
            os.path.dirname(self.bam_path),
            os.path.dirname(self.fasta_path),
            self.vcf_dir,
            self.MODEL,
        ]

    def _run(self, bam_path: str, vcf_path: str, sample_name: str = None) -> None:

        # Prepare all directories and path names for mounting
        self._mount_dirs(bam_path, vcf_path)

        # Build command
        cmd = f"singularity exec"
        cmd += f" {' '.join([f'-B {d}' for d in self.dirs])}"
        cmd += f" {self.SIF_PATH} /opt/bin/run_clair3.sh"
        cmd += f" --bam_fn={self.bam_path}"
        cmd += f" --ref_fn={self.fasta_path}"
        cmd += f" --threads='4'"
        cmd += " --platform='ont'"
        # cmd += f" --model_path=/opt/models/{self.MODEL}"
        cmd += f" --model_path={self.MODEL}"
        cmd += f" --output {self.vcf_dir}"
        cmd += " --include_all_ctgs"
        cmd += " --enable_phasing"

        # NB: We are returning all variant calls, including 0/0
        cmd += " --print_ref_calls"

        # Send all variants to the full alignment model,
        # at least theoretically maximising performance
        cmd += " --var_pct_full=1.0"
        cmd += " --ref_pct_full=1.0"

        if sample_name is not None:
            cmd += f" --sample_name={sample_name}"

        # Run
        subprocess.run(cmd, check=True, shell=True)

        # Move final VCF file produced by Clair3 to `vcf_path`
        clair3_vcf = f"{self.vcf_dir}/phased_merge_output.vcf.gz"
        shutil.copyfile(clair3_vcf, self.vcf_path)
        shutil.copyfile(clair3_vcf + ".tbi", self.vcf_path + ".tbi")


# ================================================================
# Define collection of available callers
#
# ================================================================


# Note, they are already initialised
caller_collection = {"bcftools": BcfTools, "clair3_singularity": Clair3Singularity}
