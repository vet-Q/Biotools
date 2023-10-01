import os
import click
from itertools import product
from nomadic.pipeline.cli import experiment_options, barcode_option
from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.references import PlasmodiumFalciparum3D7
from nomadic.lib.process_vcfs import bcftools_concat, bcftools_sort, bcftools_index
from nomadic.pipeline.calling.callers import caller_collection
from .downsample import BamDownSampler


@click.command(short_help="Call variants across targets.")
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
    "--downsample", is_flag=True, help="Downsample BAM before calling."
)
@click.option(
    "-r",
    "--reads",
    type=int,
    default=None,
    multiple=True,  # should allow setting multiple read depths by repeating the flag
    help="If --downsample invoked, number of reads to downsample to. Can be passed multiple times, e.g. -r 10 -r 50 -r 100.",
)
@click.option(
    "-i",
    "--iterations",
    type=int,
    default=5,
    help="If --downsample invoked, number of times to iterate downsampling for each number of reads.",
)
def call(expt_dir, config, barcode, method, downsample, reads, iterations):
    """
    Call variants, optionally with downsampling.

    """

    # Downsampling
    if downsample:
        # Check a number of reads have bee passed
        if not reads:
            raise ValueError("If --dowsample invoked, musts pass interger to -r.")

        # Run
        call_with_downsample(expt_dir, config, barcode, method, reads, iterations)

    # No dowampling
    else:
        call_all_reads(expt_dir, config, barcode, method)


def call_all_reads(expt_dir, config, barcode, method):
    """Call variants across all reads"""
    # PARSE INPUTS
    script_descrip = "NOMADIC: Call variants without downsampling"
    t0 = print_header(script_descrip)
    script_dir = f"calling/{method}"  # this will be easier to iterate overs
    params = build_parameter_dict(expt_dir, config, barcode)

    # Define reference genome
    reference = PlasmodiumFalciparum3D7()

    # Focus on a single barcode, if specified
    if "focus_barcode" in params:
        params["barcodes"] = [params["focus_barcode"]]

    # ITERATE over barcodes
    for barcode in params["barcodes"]:

        # Define input and output directory
        print("-" * 80)
        print(f"Running {method} for: {barcode}")
        barcode_dir = f"{params['barcodes_dir']}/{barcode}"
        input_dir = f"{barcode_dir}/target-extraction"
        output_dir = produce_dir(barcode_dir, script_dir)

        # ITERATE over targets
        target_vcfs = []
        for target_gene in params["target_names"]:

            # Define input bam and output vcf
            bam_path = f"{input_dir}/reads.target.{target_gene}.bam"
            vcf_path = f"{output_dir}/reads.target.{target_gene}.vcf"

            # Select variant calling method
            print("Calling variants...")
            print(f"  Input: {bam_path}")
            print(f"  Ouput: {vcf_path}")
            caller = caller_collection[method]
            caller.set_files(
                bam_path=bam_path,
                vcf_path=vcf_path,
            )
            caller.set_arguments(fasta_path=reference.fasta_path)
            caller.call_variants(sample_name=barcode)  # return True / False if it worked
            target_vcfs.append(vcf_path)
            print("Done.")
            print("")

        # Concatenate VCFs for all targets
        print("Concatenating VCFs for all targets...")
        concat_vcf = f"{output_dir}/reads.all_targets.vcf.gz"
        bcftools_concat(input_vcfs=target_vcfs, O="z", output_vcf=concat_vcf)
        sorted_vcf = concat_vcf.replace(".vcf.gz", ".sorted.vcf.gz")
        bcftools_sort(input_vcf=concat_vcf, output_vcf=sorted_vcf, O="z")
        bcftools_index(sorted_vcf)
        os.remove(concat_vcf)
        print(f"  Concatenated VCF: {sorted_vcf}")
        print("Done.")
        print("")

    print_footer(t0)


def call_with_downsample(expt_dir, config, barcode, method, reads, iterations):
    """Call variants across all reads with downsampling"""
    # PARSE INPUTS
    script_descrip = "NOMADIC: Call variants without downsampling"
    t0 = print_header(script_descrip)
    script_dir = produce_dir(
        "calling", method, "downsample"
    )  # this will be easier to iterate overs
    params = build_parameter_dict(expt_dir, config, barcode)

    # Parse reads
    print(f"Number of reads to downsample to: {reads}")

    # Define reference genome
    reference = PlasmodiumFalciparum3D7()

    # Focus on a single barcode, if specified
    if "focus_barcode" in params:
        params["barcodes"] = [params["focus_barcode"]]

    # ITERATE over barcodes
    for barcode in params["barcodes"]:
        # Define input and output directory
        print("." * 80)
        print(f"Running {method} for: {barcode}")
        print("." * 80)
        barcode_dir = f"{params['barcodes_dir']}/{barcode}"
        input_dir = f"{barcode_dir}/target-extraction"
        output_dir = produce_dir(barcode_dir, script_dir)

        # Iterate over number of reads and replicates
        for n_reads, ix in product(reads, list(range(iterations))):

            # We are creating an artifical sample from this barcode
            # by downsampling
            sample_name = f"{barcode}.n{n_reads:04d}.r{ix:03d}"
            print(f"Sample name: {sample_name}")

            # Iterate over targets
            target_vcfs = []
            for target_gene in params["target_names"]:
                print(f"Target: {target_gene}")
                # Define input bam
                bam_path = f"{input_dir}/reads.target.{target_gene}.bam"

                # Instantiate downsampler; ensure enough reads
                bam_manager = BamDownSampler(bam_path)
                if bam_manager.n_reads_total < n_reads:
                    print(f"Not enough reads for {target_gene} to downsample.")
                    print(f"  No. reads total: {bam_manager.n_reads_total}")
                    print(f"  No. downsampling: {n_reads}")
                    continue

                # Create dowsample
                print("Downsampling...")
                bam_manager.create_downsample(n_reads)
                print("Done.")
                print("")

                # Define output vcf
                vcf_fn = f"{sample_name}.{target_gene}.vcf"
                vcf_path = f"{output_dir}/{vcf_fn}"

                # Select variant calling method
                print("Calling variants...")
                print(f"  Input: {bam_path}")
                print(f"  Ouput: {vcf_path}")
                caller = caller_collection[method]
                caller.set_files(
                    bam_path=bam_manager.downsampled_bam_path,  # here we pass downsampled .bam
                    vcf_path=vcf_path,
                )
                caller.set_arguments(fasta_path=reference.fasta_path)
                status = caller.call_variants(sample_name=sample_name)
                print("Done.")
                print("")

                # Remove downsampled bam
                bam_manager.remove_downsample()

                # Store target VCF, if the file exists
                # Specifically handling issues around Clair3
                if status == 1:
                    print(f"WARNING! Clair3 DID NOT GENERATE A VCF! Not appending {vcf_path} to target VCF list.")
                    continue

                target_vcfs.append(vcf_path)

            # Concatenate for `n_reads` and `ix`
            print("Concatenating VCFs for all targets...")
            concat_vcf = f"{output_dir}/{sample_name}.all_targets.vcf.gz"
            bcftools_concat(input_vcfs=target_vcfs, O="z", output_vcf=concat_vcf)
            sorted_vcf = concat_vcf.replace(".vcf.gz", ".sorted.vcf.gz")
            bcftools_sort(input_vcf=concat_vcf, output_vcf=sorted_vcf, O="z")
            bcftools_index(sorted_vcf)
            os.remove(concat_vcf)
            print(f"  Concatenated VCF: {sorted_vcf}")
            print("Done.")
            print("")

    print_footer(t0)
