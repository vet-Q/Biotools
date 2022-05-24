import click
from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.process_vcfs import bcftools_index, bcftools_concat, bcftools_merge
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
        bcftools_index(concat_vcf)

        # Store
        vcfs_to_merge.append(concat_vcf)
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
