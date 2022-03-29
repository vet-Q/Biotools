"""
NOMADIC: Functions for mapping fastq files
2021/07/12, J.Hendry

"""


import os
import sys
import subprocess
import numpy as np
import pandas as pd


def run_minimap2(fastqs, ref_fasta, output_bam, flags=""):
    """
    Pipeline to run minimap2, compress results
    to bam file, and sort.
    
    params
        fastqs : list, str, shape (n_fastqs, )
            A list of paths to .fastq files that
            will be mapped.
        ref_fasta : str
            Path to reference genome .fasta file.
        flags : str [optional]
            Additional flags to pass to minimap2.
        output_bam : str
            Path to output BAM file.

    returns
        None
        
    """
    
    cmd = "minimap2 -ax map-ont %s %s %s | " % (ref_fasta, " ".join(fastqs), flags)
    cmd += "samtools view -S -b - | "
    cmd += "samtools sort -o %s" % output_bam
    
    subprocess.run(cmd, shell=True, check=True)
    
    return None


class Mapper:
    """
    Map a set of reads in `.fastq` format to a particular
    reference genome
    
    """
    
    def __init__(self, fastq_dir, reference):
        self.fastq_dir = fastq_dir
        self.fastqs = [f"{fastq_dir}/{fn}"
                       for fn in os.listdir(fastq_dir)
                       if fn.endswith(".fastq") or fn.endswith(".fastq.gz")]
        self.ref = reference
        
    def run_minimap2(self, output_bam, flags="--eqx --MD"):
        """
        Pipeline to run minimap2, compress results
        to bam file, and sort.

        params
            fastqs : list, str, shape (n_fastqs, )
                A list of paths to .fastq files that
                will be mapped.
            ref_fasta : str
                Path to reference genome .fasta file.
            flags : str [optional]
                Additional flags to pass to minimap2.
            output_bam : str
                Path to output BAM file.

        returns
            None

        """
        
        cmd = "minimap2"
        cmd += f" -ax map-ont {self.ref.fasta_path} {' '.join(self.fastqs)}"
        cmd += f" {flags} |"
        cmd += " samtools view -S -b - |"
        cmd += " samtools sort -o %s" % output_bam
        
        subprocess.run(cmd, shell=True, check=True)

        return None
    
