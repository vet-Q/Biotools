import click
import subprocess
from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.process_vcfs import bcftools_view, bcftools_index, bcftools_merge
from nomadic.pipeline.cli import experiment_options
from nomadic.pipeline.calling.callers import caller_collection


# TODO:
# - Most likely want to *restrict* to a particular set of barcodes
# - Also need to be able to specify which truthset to merge with
# - With all these specificatiosn, we might want a truthset-specific configuration file

# ISSUES
# - Will probably need to index and maybe also bgzip vcf files in order to run bcftools merge
# - Probably indexing will be necessary to support the merge algorithm


@click.command(short_help="Call variants across targets.")
@experiment_options
@click.option(
    "-m",
    "--method",
    type=click.Choice(caller_collection),
    required=True,
    help="Variant calling method to use.",
)
def merge(expt_dir, config, method):
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
    for barcode in params["barcodes"]:
        print(f"  Barcode: {barcode}")

        # Define directory containing VCF files
        vcf_dir = f"{params['barcodes_dir']}/{barcode}/calling/{method}"

        vcfs_to_merge = []
        print("  Compressing and indexing...")
        for target_name in params["target_names"]:

            # Compress the VCF (necessary for merging)
            input_vcf = f"{vcf_dir}/reads.target.{target_name}.vcf"
            output_vcf = f"{input_vcf}.gz"
            bcftools_view(input_vcf=input_vcf, O="z", output_vcf=output_vcf)

            # Index the VCF (also necessary for merging)
            bcftools_index(output_vcf)

            # Store
            vcfs_to_merge.append(output_vcf)
        print("Done.")
        print("")

    # Merge across all barcodes and targets
    print("Merging across all barcodes and targets...")
    print(f"  Number of .vcf files to merge: {len(vcfs_to_merge)}")
    merged_vcf = f"{output_dir}/reads.target.merged.{method}.vcf"
    bcftools_merge(input_vcfs=vcfs_to_merge, output_vcf=merged_vcf)
    print(f"  Writing to: {merged_vcf}")
    print("Done.")
    print("")
    print_footer(t0)
