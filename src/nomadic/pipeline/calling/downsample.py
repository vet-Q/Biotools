import os
import uuid
import random
import subprocess


class BamDownSampler:
    """
    Class to manage the downsampling of a `.bam` file

    """

    def __init__(self, bam_path):
        # Sanity check
        if not os.path.isfile(bam_path):
            raise FileNotFoundError(f"No .bam file found at: {bam_path}.")

        # Parse .bam location
        self.bam_path = bam_path
        self.bam_dir = os.path.dirname(bam_path)

        # Compute attributes
        self.n_reads_total = self._get_n_reads_total(bam_path)

        # Store downsampled bam information
        self.downsampled_bam_path = None

    def _get_n_reads_total(self, bam_path):
        """Compute total number of reads in a bam file"""

        cmd = f"samtools view -c {bam_path}"
        result = subprocess.run(cmd, shell=True, capture_output=True)
        n_total = int(result.stdout.decode("utf-8"))

        return n_total

    def _gen_unique_name(self):
        """Generate a unique name for the downsampled .bam file"""

        fn_id = uuid.uuid4().hex

        return f"{self.bam_dir}/temp.{fn_id}.bam"

    def create_downsample(self, n_reads, verbose=True):
        """Downsample a .bam to a specific number of reads (approximately)"""

        # Prepare an output file name
        self.downsampled_bam_path = self._gen_unique_name()

        # Compute fraction to sample
        frac = n_reads / self.n_reads_total
        if verbose:
            print(f"Downsampling to {n_reads} reads.")
            print(f"No. reads total: {self.n_reads_total}")
            print(f"Fraction to sample: {frac}")

        if frac <= 0:
            raise ValueError("Fraction of .bam to sample must be greater than zero.")

        # Prepare sub-sample tag, which requires a random seed
        rint = random.randint(0, 10 ** 9)
        stag = f"{rint}{str(frac)[1:]}"

        # Downsample
        cmd = f"samtools view"
        cmd += f" -s {stag} {self.bam_path}"
        cmd += f" -o {self.downsampled_bam_path}"
        if verbose:
            print(f"Running downsampling: {cmd}")
        subprocess.run(cmd, shell=True)

        # Index
        cmd = f"samtools index {self.downsampled_bam_path}"
        subprocess.run(cmd, shell=True)

        return None

    def remove_downsample(self):
        """Remove the most recently created downsample and it's index"""

        os.remove(self.downsampled_bam_path)
        os.remove(f"{self.downsampled_bam_path}.bai")

        return None