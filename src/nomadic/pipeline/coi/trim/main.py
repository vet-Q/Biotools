# Create a pairwise identity matrix for reads
# overlapping a target gene
# 2022/12/07, JHendry


import os
import pysam
import subprocess
import uuid
import numpy as np

from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.generic import produce_dir, print_header, print_footer
from nomadic.lib.process_bams import samtools_index, bedtools_intersect
from nomadic.lib.process_fastqs import load_fastq_read_info
from .targets import TARGET_COLLECTION


# --------------------------------------------------------------------------------
# (1) Filter BAM file
#
# --------------------------------------------------------------------------------


def bedtools_intersect_with_region(input_a, chrom, start, end, args, output):
    """
    Use `bedtools` to interesct an input file `input_a` with a
    region defined by a chromosome, start and end position

    """
    
    # Create a temporary BED file defining the region
    tmp_bed = f"tmp.{str(uuid.uuid4())[:10]}.bed"
    with open(tmp_bed, "w") as bed:
        bed.write(f"{chrom}\t{start}\t{end}\n")

    # Itersect
    bedtools_intersect(
        input_a=input_a,
        input_b=tmp_bed,
        args=args,
        output=output
    )

    # Clean
    os.remove(tmp_bed)
    

def filter_bam(
    input_bam,
    filtered_bam,
    min_read_length=2000,
    max_read_length=4000,
    min_read_qual=20):
    """
    Filter a bam file based on read length

    TODO:
    - I want to ensure these are only *primary* alignments

    """

    e =  f"length(seq) > {min_read_length}"
    e += f" && length(seq) < {max_read_length}"
    e += f" && avg(qual) > {min_read_qual}"  # could filter by read quality

    cmd = "samtools view -h"
    cmd += f" -e '{e}'"
    cmd += f" {input_bam}"
    cmd += f" -o {filtered_bam}"  # I think

    subprocess.run(cmd, check=True, shell=True)


# --------------------------------------------------------------------------------
# (2) BAM clipping and conversion to FASTQ
#
# --------------------------------------------------------------------------------


def clip_alignment(alignment, chrom, start, end):
    """
    Return sequence and quality score of of portion of
    alignment that falls within region defined by
    `chrom`, `start` and `end`
    
    params
        alignment: pysam.AlignedSegment
        start: int
        end: int
        chrom: str
        
    returns
        query_seq: str, trimmed sequence of alignment query.
        query_qual: str, trimmed quality of alignment query.
    
    """
    
    # Prepare outputs
    query_seq = ""
    query_qual = ""
    
    # Check if right chromosome
    if not alignment.reference_name == chrom:
        return query_seq, query_qual
    
    # Annoyingly, must convert to qualities to ASCII
    PYSAM_QSEQ = alignment.query_sequence
    PYSAM_QQUAL = "".join([chr(v) for v in (np.array(alignment.query_qualities) + 33)])
    assert len(PYSAM_QSEQ) == len(PYSAM_QQUAL)
        
    initialised = False
    for query_index, ref_index in alignment.aligned_pairs:

        # Skip any soft clipping
        if not initialised:
            if ref_index is None:
                continue
            else:
                initialised = True
                inside = start <= ref_index

        # Skip if it is a deletion
        if query_index is None:
            continue

        # If insertion, store if you are inside the region
        if ref_index is None:
            if not inside:
                continue
            query_seq += PYSAM_QSEQ[query_index]
            query_qual += PYSAM_QQUAL[query_index]
            continue

        if ref_index >= end:  # you have left
            break

        if start <= ref_index:
            query_seq += PYSAM_QSEQ[query_index]
            query_qual += PYSAM_QQUAL[query_index]
            inside = True
            
    assert len(query_seq) == len(query_qual)
    
    return query_seq, query_qual


def clip_bam_to_fastq(
    input_bam,
    output_fastq,
    chrom, 
    start, 
    end):
    """
    Clip alignments in a BAM file based on a region defined
    by a chromosome, start and end position
    
    """

    with open(output_fastq, "w") as fastq:
        with pysam.AlignmentFile(input_bam, "r") as bam:
            for alignment in bam:

                # Extract clipped sequence and qualities
                clipped_seq, clipped_qual = clip_alignment(alignment, chrom, start, end)

                # Write
                record =  f"@{alignment.query_name}\n"
                record += f"{clipped_seq}\n"
                record += f"+\n"
                record += f"{clipped_qual}\n"
                fastq.write(record)




# --------------------------------------------------------------------------------
# Main script
#
# --------------------------------------------------------------------------------


def main(expt_dir, config, barcode, target_gene):
    """
    Filter and trim all reads in a BAM file to overlap a `target_gene` 
    and span `start` and `end` positions; then convert to FASTQ
    
    """

    # PARSE INPUTS
    script_descrip = "NOMADIC: Trim target BAM file for COI analysis"
    t0 = print_header(script_descrip)
    script_dir = "coi"
    params = build_parameter_dict(expt_dir, config, barcode)

    target = TARGET_COLLECTION[target_gene]
    print("User inputs:")
    print(f"  Target: {target.name}")
    print(f"  Chrom: {target.chrom}")
    print(f"  Start: {target.start}")
    print(f"  End: {target.end}")
    print("Done.\n")

    # Focus on a single barcode, if specified
    if "focus_barcode" in params:
        params["barcodes"] = [params["focus_barcode"]]

    # ITERATE
    print("Iterating over barcodes...")
    for barcode in params["barcodes"]:
        print("." * 80)
        print(f"Barcode: {barcode}")
        print("." * 80)

        # DIRECTORIES
        # Input BAM
        barcode_dir = f"{params['barcodes_dir']}/{barcode}"
        bam_path = f"{barcode_dir}/target-extraction/reads.target.{target_gene}.bam"
        # Spanning [start, end] BAM
        coi_dir = produce_dir(barcode_dir, script_dir)
        bam_complete_path = f"{coi_dir}/reads.target.{target_gene}.complete.bam"

        print("Restricting to reads overlapping target region...")
        bedtools_intersect_with_region(
            input_a=bam_path,
            chrom=target.chrom, 
            start=target.start, 
            end=target.end,
            args="-F 1.0",
            output=bam_complete_path
        )
        print("Done.\n")

        # Quality and read length filtered BAM
        bam_filtered_path = bam_complete_path.replace(".bam", ".filtered.bam")

        print("Filtering BAM by read length and mean quality...")
        filter_bam(
            input_bam=bam_complete_path,
            filtered_bam=bam_filtered_path,
            min_read_length=2000,
            max_read_length=4000
        )
        samtools_index(bam_filtered_path)
        print("Done.\n")

        # FASTQ
        fastq_dir = produce_dir(coi_dir, "fastq_clipped")
        fastq_path = f"{fastq_dir}/reads.target.{target_gene}.clipped.fastq"

        print("Clipping and converting to FASTQ...")
        clip_bam_to_fastq(
            input_bam=bam_filtered_path,
            output_fastq=fastq_path,
            chrom=target.chrom,
            start=target.start,
            end=target.end
        )
        print("Done.\n")

        print("Writing read information summary CSV...")
        read_df = load_fastq_read_info(fastq_path)
        read_df.to_csv(fastq_path.replace(".fastq", ".csv"), index=False)
        print("Done.\n")
    print("Done.\n")

    print_footer(t0)
