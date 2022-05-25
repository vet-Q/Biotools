import os
import click
import numpy as np
from itertools import product
from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.process_vcfs import bcftools_sort, bcftools_index, bcftools_concat, bcftools_merge
from nomadic.pipeline.cli import experiment_options
from nomadic.pipeline.calling.callers import caller_collection


# TODO:
# - Most likely want to *restrict* to a particular set of barcodes
# - Also need to be able to specify which truthset to merge with
# - Also need to make sure edges are aligned, or else will include variants impossible to call / find truth of
# - With all these specification, we might want a truthset-specific configuration file

# ISSUES
# - Will probably need to index and maybe also bgzip vcf files in order to run bcftools merge [RESOLVED]
# - Probably indexing will be necessary to support the merge algorithm
# - How to set the names of VCF files?
#   - Best would be that all the VCF files called for a particular barcode are given a sample
#     name that corresponds to that barcode


# ================================================================
# Main script interface
#
# ================================================================


@click.command(short_help="Call variants across targets.")
@experiment_options
@click.option(
    "-m",
    "--method",
    type=click.Choice(caller_collection),
    required=True,
    help="Variant calling method to use.",
)
@click.option(
    "--downsample", 
    is_flag=True, 
    help="Produce an overview across all barcodes."
)
@click.option(
    "-r",
    "--reads",
    type=int,
    default=None,
    multiple=True, # should allow setting multiple read depths by repeating the flag
    help="If --downsample invoked, number of reads to downsample to. Can be passed multiple times, e.g. -r 10 -r 50 -r 100."
)
@click.option(
    "-i",
    "--iterations",
    type=int,
    default=10,
    help="If --downsample invoked, number of times to iterate downsampling for each number of reads."
)
def merge(expt_dir, config, method, downsample, reads, iterations):
    """ 
    Merge variant calls (VCFs) across all barcodes for a specific method,
    also merging in truth set data

    Optionally merge VCFs produced through downsampling
    
    """
    if downsample:
        merge_with_downsample(expt_dir, config, method, reads, iterations)
    else:
        merge_without_downsample(expt_dir, config, method)
    

# ================================================================
# Merge, not looking at downsampled runs
#
# ================================================================


def merge_without_downsample(expt_dir, config, method):
    """Call variants across all reads"""
    # PARSE INPUTS
    script_descrip = (
        "TRUTHSET: Merge all variants across an experiment, for a given method"
    )
    t0 = print_header(script_descrip)
    script_dir = f"calling/{method}"  # this will be easier to iterate overs
    params = build_parameter_dict(expt_dir, config, barcode=None)

    # Defie output directory
    output_dir = produce_dir(params["nomadic_dir"], script_dir)

    # Will probably want to limit to specific barcoes at some point

    # ITERATE over barcodes
    print("Preparing to merge...")
    vcfs_to_merge = []
    for barcode in params["barcodes"]:
        print(f"  Barcode: {barcode}")

        # Define directory containing VCF files
        vcf_dir = f"{params['barcodes_dir']}/{barcode}/calling/{method}"

        # Concatenating
        print("  Concatenating across targets...")
        vcfs_to_concat = [
            f"{vcf_dir}/reads.target.{target_name}.vcf"
            for target_name in params["target_names"]
        ]
        concat_vcf = f"{vcf_dir}/reads.target.concatenated.vcf.gz"
        bcftools_concat(input_vcfs=vcfs_to_concat, O="z", output_vcf=concat_vcf)
        sorted_vcf = concat_vcf.replace(".vcf.gz", ".sorted.vcf.gz")
        bcftools_sort(input_vcf=concat_vcf, output_vcf=sorted_vcf, O="z")
        bcftools_index(sorted_vcf)

        # Store
        vcfs_to_merge.append(sorted_vcf)
        print("Done.")
        print("")

    # Merge across all barcodes and targets
    print("Merging across all barcodes and targets...")
    print(f"  Number of .vcf files to merge: {len(vcfs_to_merge)}")
    merged_vcf = f"{output_dir}/reads.target.merged.{method}.vcf.gz"
    bcftools_merge(input_vcfs=vcfs_to_merge, output_vcf=merged_vcf, O="z")
    bcftools_index(merged_vcf)
    print(f"  Writing to: {merged_vcf}")
    print("Done.")
    print("")

    # MERGE WITH THE TRUTHSET
    print("Merging with the truthset...")
    METHOD = "mafft"
    truthset_vcf = f"resources/truthsets/{METHOD}/vcfs/target_genes.concat.vcf.gz"
    final_vcf = merged_vcf.replace(f"merged.{method}", f"merged.{method}.truthset.{METHOD}")
    bcftools_merge(
        input_vcfs=[merged_vcf, truthset_vcf], 
        output_vcf=final_vcf
    )
    print(f" Written to: {final_vcf}")
    print("Done.")
    print("")

    print_footer(t0)


# ================================================================
# Merge downsampled runs
#
# ================================================================


def merge_with_downsample(expt_dir, config, method, reads, iterations):
    """Call variants across all reads"""
    # PARSE INPUTS
    script_descrip = (
        "TRUTHSET: Merge all variants across an experiment, for a given method"
    )
    t0 = print_header(script_descrip)
    script_dir = f"calling/{method}"  # this will be easier to iterate overs
    params = build_parameter_dict(expt_dir, config, barcode=None)

    # Iteration specific parsing
    n_iterations = iterations
    iterations = range(n_iterations)
    print("Merging downsampled VCFs across:")
    print(f"  No. barcodes: {len(params['barcodes'])}")
    print(f"  Genes: {', '.join(params['target_names'])}")
    print(f"  Downsampled to X reads: {', '.join([str(r) for r in reads])}")
    print(f"  No. replicates: {n_iterations}")
    print("")

    # Defie output directory
    output_dir = produce_dir(params["nomadic_dir"], script_dir)

    # Will probably want to limit to specific barcoes at some point

    # ITERATE over barcodes
    print("Preparing to merge...")
    vcfs_to_merge = []
    for barcode in params["barcodes"]:
        print(f"  Barcode: {barcode}")

        # Define directory containing VCF files
        vcf_dir = f"{params['barcodes_dir']}/{barcode}/calling/{method}/downsample"
        
        # Iterate over 'experiments' = barcode/replicate/iteration and concatenate
        for n_reads, ix in product(reads, iterations):

            # Get all genes for that experiment
            vcfs_to_concat = [
                f"{vcf_dir}/{barcode}.n{n_reads:04d}.r{ix:03d}.{gene_name}.vcf"
                for gene_name in params["target_names"]
            ]

            # Concatenate
            concat_vcf = f"{vcf_dir}/reads.target.n{n_reads:04d}.r{ix:03d}.concatenated.vcf.gz"
            bcftools_concat(input_vcfs=vcfs_to_concat, O="z", output_vcf=concat_vcf)
            sorted_vcf = concat_vcf.replace(".vcf.gz", ".sorted.vcf.gz")
            bcftools_sort(input_vcf=concat_vcf, output_vcf=sorted_vcf, O="z")
            bcftools_index(sorted_vcf)

            # Store for merging across all barcodes
            vcfs_to_merge.append(sorted_vcf)
    
    print("\n".join(vcfs_to_merge))

    # Merge across all barcodes and targets
    # REPETITIVE BELOW
    print("Merging across all barcodes and targets...")
    print(f"  Number of .vcf files to merge: {len(vcfs_to_merge)}")
    merged_vcf = f"{output_dir}/reads.target.downsample.merged.{method}.vcf.gz"
    bcftools_merge(input_vcfs=vcfs_to_merge, output_vcf=merged_vcf, O="z")
    bcftools_index(merged_vcf)
    print(f"  Writing to: {merged_vcf}")
    print("Done.")
    print("")

    # MERGE WITH THE TRUTHSET
    print("Merging with the truthset...")
    METHOD = "mafft"
    truthset_vcf = f"resources/truthsets/{METHOD}/vcfs/target_genes.concat.vcf.gz"
    final_vcf = merged_vcf.replace(f"merged.{method}", f"merged.{method}.truthset.{METHOD}")
    bcftools_merge(
        input_vcfs=[merged_vcf, truthset_vcf], 
        output_vcf=final_vcf
    )
    print(f" Written to: {final_vcf}")
    print("Done.")
    print("")

    print_footer(t0)