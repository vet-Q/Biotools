import os
import subprocess
import docker
import click

from nomadic.lib.generic import produce_dir, print_header, print_footer
from nomadic.lib.parsing import build_parameter_dict
from nomadic.pipeline.cli import experiment_options, barcode_option
from nomadic.pipeline.calling.callers import caller_collection
from nomadic.lib.references import PlasmodiumFalciparum3D7


REFERENCE = PlasmodiumFalciparum3D7()
QUERY_SUBSTRING = "all_targets"


# ================================================================
# hap.py interface
#
# ================================================================


class HappyByDocker:

    IMAGE = "jmcdani20/hap.py:v0.3.12"

    def __init__(self):
        """
        Build command and run hap.py via Docker

        """

        # Set the client
        self.client = docker.from_env()

    def set_arguments(
        self,
        truth_vcf_path,
        query_vcf_path,
        reference_path,
        bed_path,
        output_dir,
        threads=6,
        happy_prefix="happy.out",
    ):
        """
        Set arguments required to run hap.py

        Paths are partitioned into directory and file name to allow
        for mounting

        """

        # TRUTH VCF
        self.truth_vcf_path = os.path.abspath(truth_vcf_path)
        self.truth_vcf_dir = os.path.dirname(self.truth_vcf_path)
        self.truth_vcf_fn = os.path.basename(self.truth_vcf_path)

        # QUERY VCF
        self.query_vcf_path = os.path.abspath(query_vcf_path)
        self.query_vcf_dir = os.path.dirname(self.query_vcf_path)
        self.query_vcf_fn = os.path.basename(self.query_vcf_path)

        # BED FILE
        self.bed_path = os.path.abspath(bed_path)
        self.bed_dir = os.path.dirname(self.bed_path)
        self.bed_fn = os.path.basename(self.bed_path)

        # REFERENCE FASTA
        self.ref_path = os.path.abspath(reference_path)
        self.ref_dir = os.path.dirname(self.ref_path)
        self.ref_fn = os.path.basename(self.ref_path)

        # OUTPUT DIRECTORY
        self.output_dir = os.path.abspath(output_dir)
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        # OUTPUT PREFIX
        self.output_prefix = f"{happy_prefix}"

        # THREADS
        self.threads = threads

        # Collect directories
        self.dirs = [
            self.truth_vcf_dir,
            self.query_vcf_dir,
            self.bed_dir,
            self.ref_dir,
            self.output_dir,
        ]

    def run(self, stratification=None):
        """
        Run hap.py via Docker

        """

        # Define volumes
        self.volumes = [f"{v}:{v}" for v in self.dirs]

        # Define command
        cmd = "/opt/hap.py/bin/hap.py"
        cmd += f" {self.truth_vcf_path}"
        cmd += f" {self.query_vcf_path}"
        cmd += f" -r {self.ref_path}"
        cmd += f" -f {self.bed_path}"
        cmd += f" -o {self.output_dir}/{self.output_prefix}"
        cmd += " --engine=vcfeval"
        cmd += f" --threads {self.threads}"
        if stratification is not None:
            stratification_path = os.path.abspath(stratification)
            cmd += f"  --stratification {stratification_path}"

            stratification_dir = os.path.dirname(stratification_path)
            if stratification_dir not in self.dirs:
                self.volumes.append(f"{stratification_dir}:{stratification_dir}")

        # Run
        output = self.client.containers.run(
            image=self.IMAGE, command=cmd, volumes=self.volumes
        )

        return output


class HappyBySingularity:

    # Prepared for Raven
    SIF_PATH = "/u/jash/containers/hap.py_v0.3.12.sif"
    SIF_CMD = "/opt/hap.py/bin/hap.py"
    BIND_DIRS = True

    def __init__(self):
        pass

    def set_arguments(
        self,
        truth_vcf_path,
        query_vcf_path,
        reference_path,
        bed_path,
        output_dir,
        threads=6,
        happy_prefix="happy.out",
    ):
        """
        Set arguments required to run hap.py

        Paths are partitioned into directory and file name to allow
        for mounting

        """

        # TRUTH VCF
        self.truth_vcf_path = os.path.abspath(truth_vcf_path)
        self.truth_vcf_dir = os.path.dirname(self.truth_vcf_path)
        self.truth_vcf_fn = os.path.basename(self.truth_vcf_path)

        # QUERY VCF
        self.query_vcf_path = os.path.abspath(query_vcf_path)
        self.query_vcf_dir = os.path.dirname(self.query_vcf_path)
        self.query_vcf_fn = os.path.basename(self.query_vcf_path)

        # BED FILE
        self.bed_path = os.path.abspath(bed_path)
        self.bed_dir = os.path.dirname(self.bed_path)
        self.bed_fn = os.path.basename(self.bed_path)

        # REFERENCE FASTA
        self.ref_path = os.path.abspath(reference_path)
        self.ref_dir = os.path.dirname(self.ref_path)
        self.ref_fn = os.path.basename(self.ref_path)

        # OUTPUT DIRECTORY
        self.output_dir = os.path.abspath(output_dir)
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        # OUTPUT PREFIX
        self.output_prefix = f"{happy_prefix}"

        # THREADS
        self.threads = threads

        # Collect directories
        self.dirs = [
            self.truth_vcf_dir,
            self.query_vcf_dir,
            self.bed_dir,
            self.ref_dir,
            self.output_dir,
        ]

    def run(self, stratification=None):
        """
        Run hap.py via Singularity

        """

        # Add stratification to mounting
        if stratification is not None:
            stratification_path = os.path.abspath(stratification)
            stratification_dir = os.path.dirname(stratification_path)
            if stratification_dir not in self.dirs:
                self.dirs.append(stratification_dir)

        # Define command
        cmd = f"singularity exec"
        if self.BIND_DIRS:
            cmd += f" {' '.join([f'-B {d}' for d in self.dirs])}"
        cmd += f" {self.SIF_PATH} {self.SIF_CMD}"
        cmd += f" {self.truth_vcf_path}"
        cmd += f" {self.query_vcf_path}"
        cmd += f" -r {self.ref_path}"
        cmd += f" -f {self.bed_path}"
        cmd += f" -o {self.output_dir}/{self.output_prefix}"
        cmd += " --engine=vcfeval"
        cmd += f" --threads {self.threads}"
        if stratification is not None:
            cmd += f"  --stratification {stratification_path}"

        # Run
        subprocess.run(cmd, check=True, shell=True)

        return cmd


happy_callers = {"docker": HappyByDocker, "singularity": HappyBySingularity}


# ================================================================
# Main script
#
# ================================================================


@click.command(short_help="Comparing VCFs with hap.py.")
@experiment_options
@barcode_option
@click.option(
    "-m",
    "--method",
    type=click.Choice(caller_collection),
    required=True,
    help="Variant calling method to use.",
)
@click.option(
    "-t",
    "--truth_vcf",
    type=click.Path(),
    required=True,
    help="VCF file to be taken as truth.",
)
@click.option(
    "-f",
    "--bed_path",
    type=click.Path(),
    required=True,
    help="Path to BED file defining regions in which"
    " all variants in `truth_vcf` will be considered true positives.",
)
@click.option(
    "-s",
    "--stratification",
    type=click.Path(),
    required=False,
    default=None,
    help="Path to TSV file for stratifications.",
)
@click.option(
    "--downsample",
    is_flag=True,
    help="Compare against downsampled vcfs (e.g. in /downsample).",
)
@click.option(
    "-h",
    "--happy_caller",
    type=click.Choice(happy_callers),
    required=False,
    default="docker",
    help="Run hap.py by Docker or Singularity.",
)
def cfhappy(
    expt_dir,
    config,
    barcode,
    method,
    truth_vcf,
    bed_path,
    stratification,
    downsample,
    happy_caller,
):
    """
    Run a comparison with `hap.py` for a specific `barcode` and `method` from an
    `expt_dir` against a `truth_vcf`

    """
    # PARSE INPUTS
    print("Parsing input arguments...")
    script_descrip = "TRUTHSET: Compare VCFs from a given barcode with a truth set."
    t0 = print_header(script_descrip)
    params = build_parameter_dict(expt_dir, config, barcode=barcode)

    # Define critical directories
    script_dir = f"{params['focus_barcode']}/cfhappy/{method}"
    output_dir = produce_dir(params["barcodes_dir"], script_dir)
    input_dir = output_dir.replace("cfhappy", "calling")
    if downsample:
        input_dir += "/downsample"

    # Locate query VCFs
    query_vcfs = [
        f
        for f in os.listdir(input_dir)
        if QUERY_SUBSTRING in f and f.endswith(".vcf.gz")
    ]

    print("  Directories")
    print(f"    Input directory: {input_dir}")
    print(f"    Output directory: {output_dir}")
    print("  Files")
    print(f"    Truth VCF: {truth_vcf}")
    print(f"    Reference FASTA: {REFERENCE.fasta_path}")
    print(f"    BED path: {bed_path}")
    if stratification is not None:
        print(f"    Stratifications TSV: {stratification}")
    print(f"    Comparing against downsample? {downsample}")
    print(f"    Query VCF substring: {QUERY_SUBSTRING}")
    print(f"    No. VCFs found to compare: {len(query_vcfs)}")
    print("Done.")
    print("")

    # Prepare hap.py API
    print("Running hapy.py...")
    happy = happy_callers[happy_caller]()

    # Iterate over query VCFs
    for query_vcf in query_vcfs:
        print(f"  Query VCF: {query_vcf}")
        happy.set_arguments(
            truth_vcf_path=truth_vcf,
            query_vcf_path=f"{input_dir}/{query_vcf}",
            reference_path=REFERENCE.fasta_path,
            bed_path=bed_path,
            output_dir=output_dir,
            happy_prefix=query_vcf.replace(".vcf", ""),
        )
        output = happy.run(stratification=stratification)
        print(f"  Outputs written to: {output_dir}")
        print("Done.")  # TODO: better way to check outputs
        print("")

    print_footer(t0)
