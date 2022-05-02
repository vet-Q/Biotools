import click
from nomadic.pipeline.cli import experiment_options, barcode_option
from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.references import PlasmodiumFalciparum3D7
from nomadic.pipeline.calling.callers import caller_collection


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
    "--downsample", is_flag=True, help="Produce an overview across all barcodes."
)
def call(expt_dir, config, barcode, method, downsample):
    """
    Run analyses of a specific set of amplicon targets

    TODO:
    - Probably want an additional argument here specifying targets
    - Would add flexibility

    """
    # PARSE INPUTS
    script_descrip = "NOMADIC: Map .fastq files to Plasmodium falciparum"
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
        print("-"*80)
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
            caller.set_arguments(
                fasta_path=reference.fasta_path
            )
            caller.call_variants()
            print("Done.")
            print("")
    print_footer(t0)


def call_with_downsample(expt_dir, config, barcode):
    """
    Call variants for a P.f. sample

    TODO:
    - INCOMPLETE

    Steps
    - Load params
    - Run with or without downsampling;
    - Call variants across all extracted targets

    """

    # PARSE INPUTS
    script_descrip = "NOMADIC: Map .fastq files to Plasmodium falciparum"
    t0 = print_header(script_descrip)
    script_dir = "calling"
    params = build_parameter_dict(expt_dir, config, barcode)

    # Focus on a single barcode, if specified
    if "focus_barcode" in params:
        params["barcodes"] = [params["focus_barcode"]]

    # Define reference genomes
    references = [
        PlasmodiumFalciparum3D7(),
    ]

    # ITERATE over barcodes
    for barcode in params["barcodes"]:
        
        # Define input and output directory
        print("-"*80)
        print(f"Running longshot for: {barcode}")
        barcode_dir = f"{params['barcodes_dir']}/{barcode}"
        input_dir = f"{barcode_dir}/target-extraction"
        output_dir = produce_dir(barcode_dir, script_dir)

        # ITERATE over targets
        for target_gene in params["target_names"]:

            # Define input bam
            bam_path = f"{input_dir}/reads.target.{target_gene}.bam"

            # Instantiate downsampler; ensure enough reads
            bam_manager = BamDownSampler(bam_path)
            if bam_manager.n_reads_total < args.downsample_reads:
                print(f"Not enough reads for {target_gene} to downsample.")
                print(f"  No. reads total: {bam_manager.n_reads_total}")
                print(f"  No. downsampling: {args.downsample_reads}")
                continue

            # ITERATE over replicates
            for ix in range(args.iterations):

                # Create dowsample
                bam_manager.create_downsample(args.downsample_reads)

                # Define output vcf
                vcf_fn = f"{args.barcode}.{target_gene}."
                vcf_fn += f"n{args.downsample_reads:04d}.r{ix:03d}.vcf"
                vcf_path = f"{output_dir}/{vcf_fn}"

                # Select variant calling method

                # Run

                # Remove downsample



