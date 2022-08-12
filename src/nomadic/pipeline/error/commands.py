# TODO
# - Sub-directories for approach
#
# Overview
# for barcode in barcodes...
# select bam
# create .bed
# create .pileup [optionally invert]
# summarise in .csv

import click
from nomadic.pipeline.cli import experiment_options, barcode_option
from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.references import PlasmodiumFalciparum3D7
from nomadic.lib.process_gffs import load_gff, add_gff_fields
from .algorithm import ErrorAnalysisAlgorithm
from .bed import ByProteinCodingGene


@click.command(short_help="Characterise error rate.")
@experiment_options
@barcode_option
@click.option(
    "-a",
    "--approach",
    type=str, # will be from choices
    #required=True,
    help="Variant calling method to use.",
)
def error(expt_dir, config, barcode, approach):
    # PARSE INPUTS
    script_descrip = "NOMADIC: Characterise all differences from the reference genome"
    t0 = print_header(script_descrip)
    script_dir = "error"
    params = build_parameter_dict(expt_dir, config, barcode)

    # Focus on a single barcode, if specified
    if "focus_barcode" in params:
        params["barcodes"] = [params["focus_barcode"]]

    # Define reference genomes
    print("Defining reference genome")
    reference = PlasmodiumFalciparum3D7()
    print(f"  Name: {reference.name}")
    print("Done.")
    print("")

    # Load gff
    print("Loading .gff")
    print(f"  from: {reference.gff_path}")
    gff_df = load_gff(reference.gff_path)
    gff_df = add_gff_fields(gff_df)
    print(f"  No. records: {gff_df.shape[0]}")
    print("Done.")
    print("")

    # Bed builder
    bed_builder = ByProteinCodingGene(gff_df)

    # ITERATE over barcodes
    print("Iterating over barcodes...")
    for barcode in params["barcodes"]:
        print("." * 80)
        print(f"Barcode: {barcode}")
        print("." * 80)

        # Define output directory, initiate algorithm
        output_dir = produce_dir(params['barcodes_dir'], barcode, script_dir)
        algorithm = ErrorAnalysisAlgorithm(reference=reference, output_dir=output_dir)

        # ITERATE over targets
        for target_id, target_name in params["name_dt"].items():
            print(f"  Target: {target_id} = {target_name}")

            # Run algorithm to produce summary data rfames
            bam_path = f"{params['barcodes_dir']}/{barcode}/target-extraction/reads.target.{target_name}.bam"
            algorithm.set_target(target_id)
            algorithm.create_target_bed(bed_builder)
            algorithm.create_target_mpileup(bam_path=bam_path)
            mutation_df, indel_df = algorithm.get_mpileup_summary()

            # Annotate
            mutation_df.insert(0, "ID", target_id)
            mutation_df.insert(1, "gene_name", params["name_dt"][target_id])
            indel_df.insert(0, "ID", target_id)
            indel_df.insert(1, "gene_name", params["name_dt"][target_id])
            
            # Save
            mutation_df.to_csv(algorithm.pileup_path.replace(".mpileup", ".nt_error.csv"))
            indel_df.to_csv(algorithm.pileup_path.replace(".mpileup", ".indel_lengths.csv"))



        
        



     #   # Create mutation and indel data frames
#         print("    Analysing basecalls...")
#         mutation_df, indel_df 
        
#         # Annotate
#         mutation_df.insert(0, "ID", target_id)
#         mutation_df.insert(1, "gene_name", params["name_dt"][target_id])
#         indel_df.insert(0, "ID", target_id)
#         indel_df.insert(1, "gene_name", params["name_dt"][target_id])
        
#         # Save
#         mutation_df.to_csv(output_pileup.replace(".mpileup", ".nt_error.csv"))
#         indel_df.to_csv(output_pileup.replace(".mpileup", ".indel_lengths.csv"))
        
#         # Store
#         all_mutation_df.append(mutation_df)
#         all_indel_df.append(indel_df)   