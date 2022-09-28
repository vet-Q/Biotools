"""
NOMADIC: Functions for processing bam files
2021/08/05, JHendry

"""

import subprocess
import pandas as pd
import numpy as np
import re




# ================================================================================ #
# Functions for running samtools utilities 
#
# ================================================================================ #




def run_samtools(utility, args, output_path=None):
    """
    Run a samtools utility
    
    params
        utility : str
            Name of the samtools utility.
        args : str
            Arguments passed to utility
        output_path : str [optional]
            Optionally redirect output to `output_path`.
            
    returns
        None
    
    """
    
    cmd = "samtools %s %s" % (utility, args)
    if output_path is not None:
        cmd += "> %s" % output_path
    subprocess.run(cmd, shell=True, check=True)
    
    return None


def samtools_view(input_bam, args, output_bam):
    """
    Run `samtools view` on an `input_bam`
    
    params
        input_bam : str
            Path to bam file.
        args : str
            Arguments to pass to `samtools view`
        output_bam : str
            Path to output bam file.
    
    returns
        None
        
    """
    
    cmd = "samtools view %s %s -o %s" % (input_bam, args, output_bam)
    subprocess.run(cmd, check=True, shell=True)
    
    return None


def samtools_index(input_bam):
    """
    Run index BAM
    
    params
        input_bam : str
            BAM file to be indexed.
    
    returns
        None

    """
    
    cmd = "samtools index %s" % input_bam
    subprocess.run(cmd, shell=True, check=True)
    
    return None


def samtools_merge(bam_files, output_bam):
    """
    Merge a collection of BAM files `bam_files`
    and write as an output BAM `output_bam`
    
    params
        bam_files : list, str, shape (n_bams, )
            A list of paths to bam files that should
            be merged.
        output_bam : str
            Path to output BAM file.
    
    returns
        None

    """
    
    cmd = "samtools merge -f %s %s" % (output_bam, " ".join(bam_files))
    subprocess.run(cmd, shell=True, check=True)

    return None


def samtools_mpileup(input_bam, ref_fasta, target_bed, output_pileup, downsample=False):
    """
    Generate a pileup over an `input_bam` file for a
    given `target_gene` and write to `output_pileup`
    
    
    params
        input_bam : str
            BAM file which will be used to generate
            pileup.
        ref_fasta : str
            Path to reference genome in .fasta format.
        target_bed : ?
            ?
        output_pileup : str
            Path to output pileup file.
    
    returns
        None
    
    """
    
    if downsample:
        pass
    
    cmd = "samtools mpileup"
    cmd += " -B %s" % input_bam
    cmd += " -f %s" % ref_fasta
    cmd += " -l %s" % target_bed
    cmd += " -Q 0"  # Minimum quality for a base to be considered
    cmd += " -aa"  # Output absolutely all positions, including unused reference sequence
    cmd += " > %s" % output_pileup 
    
    subprocess.run(cmd, shell=True, check=True)
    
    return None


def samtools_stats(input_bam, output_stats, view_args=None):
    """
    Run `samtools stats` on a given `input_bam`,
    optionally filtering first with `samtools view`
    
    params
        input_bam : str
            Path to bam file.
        output_stats : str
            Path to write output file.
        view_args : str [optional]
            Arguments to pass to `samtools view`,
            that allow for filtering before
            `samtools stats` is run. e.g., could
            be '-F 0x904' to capture only mapped reads.
            
    returns
        None
    
    """
    
    cmd = "samtools view -b %s %s" % ("" if view_args is None else view_args, input_bam)
    cmd += " | samtools stats -"
    cmd += " > %s" % output_stats
    subprocess.run(cmd, shell=True, check=True)
    
    return None


def samtools_depth(input_bam, output_path, region_bed=None):
    """
    Run `samtools depth` on a given `input_bam`, focussing
    on regions defined by a `bed_file`

    params
        input_bam : str
            Path to bam file.
        output_path : str
            Path to write output file.
        region_bed : str
            BED file defining regions over which depth
            should be calculated. [optional]
    
    returns
        None

    """

    cmd = "samtools depth"
    cmd += " -aa" # output all positions
    cmd += f" -b {region_bed}"
    cmd += f" -o {output_path}"
    cmd += f" {input_bam}"
    subprocess.run(cmd, shell=True, check=True)

    return None


# ================================================================================ #
# Functions for running bedtools utilities 
#
# ================================================================================ #




def bedtools_intersect(input_a, input_b, args, output):
    """
    Run `bedtools intersect` on an `input_bam`
    
    params
        input_a, input_b : str
            Paths to BAM/BED/GFF &c files...
        args : str
            Arguments to pass to `bedtools intersect`
            
    returns
        None
    
    """
    
    cmd = "bedtools intersect -a %s -b %s %s > %s" % (input_a, input_b, args, output)
    subprocess.run(cmd, check=True, shell=True)
    
    return None




# ================================================================================ #
# Functions for processing bam files or bam stats files
#
# ================================================================================ #




def extract_bam_MAPQ(input_bam):
    """
    Extract the mapping quality scores of a target
    BAM file.
    
    Qualtities Q are given by:
    
    Q = -10log_{10} Pr(read is mapped wrong)
    
    params
        input_bam : str
            Path to target bam file.
    
    returns
        mapqs : ndarray, int8, shape(n, )
            Mapping qualities for each alignment.
    
    """
    
    cmd = "samtools view %s | cut -f 5" % input_bam
    result = subprocess.run(cmd, shell=True, capture_output=True)
    stdout = result.stdout.decode("utf-8")
    mapqs = np.array(stdout.split(), "int8")
    
    return mapqs


def extract_bam_FLAG(input_bam):
    """
    Extract the bitwise flag column from a target
    BAM file.
    
    The FLAG column gives qualitative information
    about the nature of the alignment, i.e. whether
    it is supplementary, secondary, on the reverse
    strand, &c.
    
    params
        input_bam : str
            Path to target bam file.
    
    returns
        flags : ndarray, int16, shape(n, )
        Bitwise flags for each alignment as
        an integer.
    
    """
    
    cmd = "samtools view %s | cut -f 2" % input_bam
    result = subprocess.run(cmd, shell=True, capture_output=True)
    stdout = result.stdout.decode("utf-8")
    flags = np.array(stdout.split(), "int16")
    
    return flags


def extract_bam_stats_RL(bam_stats):
    """
    Parse the read length `RL` section of `samtool stats` output
    
    params
        bam_stats : str
            Path to a file produced by running
            `samtool stats` on a .bam file.
    returns
        _ : DataFrame, int, shape(n, 2)
            Data frame giving the read length
            histogram.
    
    """
    
    # Extract relevant lines using grep
    cmd = "cat %s | grep ^RL | cut -f 2,3" % bam_stats

    # Run command
    result = subprocess.run(cmd, shell=True, capture_output=True)

    # Parse stdout into DataFrame
    stdout = result.stdout.decode("utf-8")
    dt = {
        "read_length" : np.array(stdout.split()[:-1:2], "int"),
        "count" : np.array(stdout.split()[1::2], "int")
    }

    return pd.DataFrame(dt)


def extract_bam_stats_COV(bam_stats):
    """
    Parse the coverage `COV` section of `samtool stats` output
    
    params
        bam_stats : str
            Path to a file produced by running
            `samtool stats` on a .bam file.
            
    returns
        _ : DataFrame, int, shape(n, 2)
            Data frame giving the coverage histogram
            histogram.
    
    """
    
    # Extract relevant lines using grep
    cmd = "cat %s | grep ^COV | cut -f 3,4" % bam_stats

    # Run command
    result = subprocess.run(cmd, shell=True, capture_output=True)

    # Parse stdout into DataFrame
    stdout = result.stdout.decode("utf-8")
    dt = {
        #"coverage_bin" : np.array(stdout.split()[:-1:2], "int"),
        "coverage" : np.array(stdout.split()[:-1:2], "int"),
        "count" : np.array(stdout.split()[1::2], "int")
    }

    return pd.DataFrame(dt)


def summarise_bam_stats(input_file, view_args=None):
    """
    Calculate statistics for an `input_bam` and return
    as a dictionary
    
    params
        input_file : str
            Path to either .bam or .stats file.
        view_args : str [optional]
            Arguments to pass to `samtools view`,
            that allow for filtering before
            `samtools stats` is run. e.g., could
            be '-F 0x904' to capture only mapped reads.
    
    returns
        stat_dt : dict
            Dictionary containing summary statistics
            calculated on `input_bam`.
    
    """
    
    # Define command
    if input_file.endswith(".bam"):
        cmd = "samtools view -b %s %s" % ("" if view_args is None else view_args,
                                          input_file)
        cmd += " | samtools stats -- "
    elif input_file.endswith(".stats"):
        cmd = "cat %s" % input_file
    else:
        raise ValueError("`input_file` must end in `.bam` or `.stats`.")
    cmd += " | grep ^SN"  # limit to summary statistics
    cmd += " | cut -f 2-"

    # Run command
    results = subprocess.run(cmd, shell=True, capture_output=True)

    # Partition summary lines
    lines = results.stdout.decode("utf-8").split("\n")

    # Extract key information into dict
    target_lines = {
        "raw total sequences:": "reads_total", 
        "reads mapped:": "reads_mapped", 
        "total length:": "bases_total", 
        "bases mapped (cigar):": "bases_mapped",
        "average length:": "mean_read_length",
        "average quality:": "mean_read_qual",
        "mismatches:": "mismatches",
        "error rate:": "error_rate"
    }

    stat_dt = {}
    for k, v in target_lines.items():
        for line in lines:
            if line.startswith(k):
                r = line.split("\t")
                value = float(r[1])
                stat_dt[v] = value

    return stat_dt




# ================================================================================ #
# Functions running samtools utilities and doing some processing
#
# ================================================================================ #




def get_mapping_statistics(bam, verbose=True):
    """
    Get mapping statistics for a given `bam` file
    using samtools
    
    Statistics are computed using the bit-wise flag (second
    field in a SAM record). Relevant bits are include...
    
    0x4  Segment is unmapped
    0x100  Secondary alignment (non-primary multiple)
    0x800  Supplementary alignment (non-primary chimeric)
    
    The flag "-f" retrieves all SAM records where the bit(s)
    is/are set. "-F" retrieves all SAM records where they are
    not set.
    
    
    Some key concepts include...
    
    Linear alignment:
    An alignment of a read to a single reference sequence
    that may include insertions, deletions, skips and clipping,
    but may not include direction changes (i.e. part forward
    and part reverse strand). A linear alignment can be 
    represented in a single SAM record.
    
    Multiple mapping: 
    The alignment of a read may be ambiguous, e.g. as a 
    result of repeats. In this case the read will map to 
    multiple locations. One of these mappings will be 
    considered primary. All others have the
    *secondary* flag set.
    
    Chimeric alignment:
    An alignment of a read that cannot be represented 
    as a linear alignment. A chimeric alignment is
    represented as a set of linear alignments that do
    not have large overlaps. Typically, one alignment
    is defined as 'representative' and all the others
    have a *supplementary* flag set.

    
    params
        bam : str
            Path to target bam file.
        verbose : bool
            Print the mapping statistics to screen.
    returns
        stat_dt : dict
            A dictionary containing information on
            a set of mapping statistics.
    
    """
    
    samtools = "samtools view -c"  # return number of records
    
    flag_dt = {
        "Alignments": "",
        "Reads": "-F 0x900",  # Do not count supplementary
        "Mapped": "-F 0x904",
        "Unmapped": "-f 0x004",
        "Multi-mapped": "-f 0x100",
        "Chimera-mapped": "-f 0x800"
    }
    
    stat_dt = {}
    for stat, flag in flag_dt.items():
        result = subprocess.check_output(" ".join([samtools, flag, bam]), shell=True)
        stat_dt[stat] = int(result)
        
    if verbose:
        print("Mapping Statistics")
        print("  BAM: %s" % bam)
        n_reads = stat_dt["Reads"]
        for stat, value in stat_dt.items():
            print("  %s: %d (%.02f%%)" % (stat, value, 100*value/n_reads))
        
    return stat_dt


def combine_mapping_statistics(mapping_statistics, rescale=10**3):
    """
    Combine mapping statistics across barcodes
    
    """
    
    df = pd.DataFrame(mapping_statistics)
    
    df = df[["Barcode", "Alignments", "Reads", 
             "Mapped", "Unmapped",
             "Multi-mapped", "Chimera-mapped"]]
    
    df.sort_values("Barcode", inplace=True)
    df.reset_index(drop=True, inplace=True)
    df.index = df["Barcode"]
    df.drop("Barcode", 1, inplace=True)
    df = df / rescale
    
    return df


# ================================================================================ #
# Functions for processing the outputs of the above functions
#
# ================================================================================ #


def convert_RL_to_ndarray(rl_df):
    """
    Convert a read length data frame to an array of read lengths
    
    params
        rl_df: DataFrame, shape (n_read_lengths, 2)
            Read length dataframe, as would be returned by 
            `extract_bam_stats_RL`.
    
    returns
        _ : ndarray, shape(n_reads, )
            Array of read lengths.
    
    """
    
    assert "read_length" in rl_df.columns, "'read_length' must be a column."
    assert "count" in rl_df.columns, "'count' must be a column."
    
    rls = []
    for i, row in rl_df.iterrows():
        rls += [row["read_length"]] * row["count"]
        
    return np.array(rls)

