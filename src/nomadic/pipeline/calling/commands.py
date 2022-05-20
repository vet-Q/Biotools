import click
from nomadic.pipeline.cli import experiment_options, barcode_option
from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.references import PlasmodiumFalciparum3D7
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
def call(expt_dir, config, barcode, method, downsample, reads, iterations):
    """
    Call variants, optionally with downsampling.
    
    """
    
    # Downsampling
    if downsample:
        # Check a number of reads have bee passed
        if reads is None:
            raise ValueError("If --dowsample invoked, musts pass interger to -r.")

        # Run
        call_with_downsample(expt_dir, config, barcode, method, reads, iterations)
        
    # No dowampling
    else:
        call_all_reads(expt_dir, config, barcode, method)


def call_all_reads(expt_dir, config, barcode, method):
    """ Call variants across all reads """
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
            caller.call_variants(sample_name=barcode)
            print("Done.")
            print("")
    print_footer(t0)


def call_with_downsample(expt_dir, config, barcode, method, reads, iterations):
    """

    """
    # PARSE INPUTS
    script_descrip = "NOMADIC: Call variants without downsampling"
    t0 = print_header(script_descrip)
    script_dir = produce_dir("calling", method, "downsample")  # this will be easier to iterate overs
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

        # ITERATE over reads
        for target_gene in params["target_names"]:
            print(f"Target: {target_gene}")
            # Define input bam
            bam_path = f"{input_dir}/reads.target.{target_gene}.bam"
            for n_reads in reads:
                # Instantiate downsampler; ensure enough reads
                bam_manager = BamDownSampler(bam_path)
                if bam_manager.n_reads_total < n_reads:
                    print(f"Not enough reads for {target_gene} to downsample.")
                    print(f"  No. reads total: {bam_manager.n_reads_total}")
                    print(f"  No. downsampling: {n_reads}")
                    continue
                
                for ix in range(iterations):

                    # Create dowsample
                    print("Downsampling...")
                    bam_manager.create_downsample(n_reads)
                    print("Done.")
                    print("")

                    # Define output vcf
                    sample_name = f"{barcode}.n{n_reads:04d}.r{ix:03d}"
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
                    caller.call_variants(sample_name=sample_name)
                    print("Done.")
                    print("")

                    # Remove downsampled bam
                    bam_manager.remove_downsample()
    print_footer(t0)











# @click.group(short_help="Call variants across targets.")
# @experiment_options
# @barcode_option
# @click.option(
#     "-m",
#     "--method",
#     type=click.Choice(caller_collection),
#     required=True,
#     help="Variant calling method to use.",
# )
# @click.option(
#     "--downsample", is_flag=True, help="Produce an overview across all barcodes."
# )
# def call(expt_dir, config, barcode, method, downsample):
#     """
#     Run analyses of a specific set of amplicon targets

#     TODO:
#     - Probably want an additional argument here specifying targets
#     - Would add flexibility

#     """
#     # PARSE INPUTS
#     script_descrip = "NOMADIC: Map .fastq files to Plasmodium falciparum"
#     t0 = print_header(script_descrip)
#     script_dir = f"calling/{method}"  # this will be easier to iterate overs
#     params = build_parameter_dict(expt_dir, config, barcode)

#     # Define reference genome
#     reference = PlasmodiumFalciparum3D7()

#     # Focus on a single barcode, if specified
#     if "focus_barcode" in params:
#         params["barcodes"] = [params["focus_barcode"]]

#     # ITERATE over barcodes
#     for barcode in params["barcodes"]:

#         # Define input and output directory
#         print("-" * 80)
#         print(f"Running {method} for: {barcode}")
#         barcode_dir = f"{params['barcodes_dir']}/{barcode}"
#         input_dir = f"{barcode_dir}/target-extraction"
#         output_dir = produce_dir(barcode_dir, script_dir)

#         # ITERATE over targets
#         for target_gene in params["target_names"]:

#             # Define input bam and output vcf
#             bam_path = f"{input_dir}/reads.target.{target_gene}.bam"
#             vcf_path = f"{output_dir}/reads.target.{target_gene}.vcf"

#             # Select variant calling method
#             print("Calling variants...")
#             print(f"  Input: {bam_path}")
#             print(f"  Ouput: {vcf_path}")
#             caller = caller_collection[method]
#             caller.set_files(
#                 bam_path=bam_path,
#                 vcf_path=vcf_path,
#             )
#             caller.set_arguments(fasta_path=reference.fasta_path)
#             caller.call_variants()
#             print("Done.")
#             print("")
#     print_footer(t0)

