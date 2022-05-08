import click
import pandas as pd
from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.references import PlasmodiumFalciparum3D7
from nomadic.pipeline.cli import experiment_options
from nomadic.pipeline.calling.callers import caller_collection
from nomadic.pipeline.find.gene import Gene
from nomadic.pipeline.find.vcf import load_vcf_using_allel
from nomadic.pipeline.find.plot import MutationPanelPlot


@click.command(short_help="Find a set of target mutations.")
@experiment_options
@click.option(
    "-m",
    "--method",
    type=click.Choice(caller_collection),
    required=True,
    help="Variant calling method to use.",
)
def find(expt_dir, config, method):
    """
    Find a specific set of mutations across target amplicons

    TODO:
    - PASS COVERAGE INFORMATION -- important piece of context
      for interpreting genotyping results

    """

    # PARSE INPUTS
    script_descrip = "NOMADIC: Find mutations in P.f. amplicon data"
    t0 = print_header(script_descrip)
    script_dir = "find"
    params = build_parameter_dict(expt_dir, config, barcode=None)

    # Create output directory
    output_dir = produce_dir(params["nomadic_dir"], script_dir)

    # Define reference genomes
    reference = PlasmodiumFalciparum3D7()

    # Create a simple resistance mutation data frame
    resistance_df = params["mutations"]
    resistance_df["resistance"] = True

    # Create gene set based on targets
    gene_set = {}
    for target_id in params["target_ids"]:
        print(f"Creating gene {target_id}.")
        gene = Gene(gff_path=reference.gff_path, fasta_path=reference.fasta_path)
        gene.set_gene_id(target_id)
        gene_set[target_id] = gene

    # ITERATE over barcodes
    overall_results = []
    for barcode in params["barcodes"]:
        print(f"Iterating over: {barcode}")

        # Define the output directory for this barcode
        barcode_output_dir = produce_dir(params["barcodes_dir"], barcode, script_dir)

        # Iterate over of target genes
        target_results = []
        for target_id, gene_name in params["name_dt"].items():

            # Load VCF (probably need a method argument)
            vcf_path = f"{params['barcodes_dir']}/{barcode}/calling/{method}/reads.target.{gene_name}.vcf"
            vcf_df = load_vcf_using_allel(vcf_path)

            # Extract gene of interets
            focus_gene = gene_set[target_id]

            # Label all mutations
            mutations = [
                focus_gene.query_for_mutation(
                    chrom=row["chrom"], pos=row["pos"], ref=row["ref"], alt=row["alt"]
                )
                for _, row in vcf_df.iterrows()
            ]

            # Insert new column
            vcf_df.insert(0, "mutation", mutations)

            # Annotate with barcode and gene
            vcf_df.insert(0, "barcode", barcode)
            vcf_df.insert(1, "target_id", target_id)
            vcf_df.insert(2, "gene_name", gene_name)

            # Store
            target_results.append(vcf_df)

        # Create dataframe for the barcode
        barcode_df = pd.concat(target_results)

        # Create a resistace only-dataframe
        barcode_resistance_df = pd.merge(
            left=resistance_df,
            right=barcode_df[["gene_name", "mutation", "gt", "qual"]],
            how="left",
            on=["gene_name", "mutation"],
        ).fillna(0)
        barcode_resistance_df.insert(0, "barcode", barcode)

        # Write
        barcode_df.to_csv(f"{barcode_output_dir}/table.mutation_summary.{method}.csv")
        barcode_resistance_df.to_csv(
            f"{barcode_output_dir}/table.mutation_summary.resistance.{method}.csv"
        )

        # Store overall results
        overall_results.append(barcode_resistance_df)

    # Combine, merge with metadata, save
    overall_df = pd.concat(overall_results)
    overall_df = pd.merge(left=params["metadata"], right=overall_df, on="barcode")
    overall_df.to_csv(f"{output_dir}/table.resistance_mutations.{method}.gt.csv")

    # Plot and save
    plotter = MutationPanelPlot(overall_df)
    plotter.plot_multiple_genes(
        output_path=f"{output_dir}/plot.resistance.{method}.gt.pdf"
    )

    print_footer(t0)
