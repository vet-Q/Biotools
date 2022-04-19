from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.references import PlasmodiumFalciparum3D7
from nomadic.pipeline.find.gene import Gene
from nomadic.pipeline.find.vcf import load_vcf_as_df


def main(expt_dir, config, barcode, method):
    """
    Find a specific set of mutations across target amplicons

    """

    # PARSE INPUTS
    script_descrip = "NOMADIC: Find mutations in P.f. amplicon data"
    t0 = print_header(script_descrip)
    script_dir = "find"
    params = build_parameter_dict(expt_dir, config, barcode)

    # Create output directory
    output_dir = produce_dir(params["nomadic_dir"], script_dir)

    # Focus on a single barcode, if specified
    if "focus_barcode" in params:
        params["barcodes"] = [params["focus_barcode"]]

    # Define reference genomes
    reference = PlasmodiumFalciparum3D7()

    # Create gene set based on targets
    # - Do this once to reduce IO, rather tha for each barcode
    gene_set = {}
    for target_id in params["target_ids"]:
        print(f"Creating gene {target_id}.")
        gene = Gene(gff_path=reference.gff_path, fasta_path=reference.fasta_path)
        gene.set_gene_id(target_id)
        gene_set[target_id] = gene

    # ITERATE over barcodes
    dfs = []
    for barcode in params["barcodes"]:
        print(f"Iterating over: {barcode}")
        
        # Go for one directory for the entire experiment (?)
        # Or save merging as a separate script (?)
        # Create output directory
        #output_dir = produce_dir(params['barcodes_dir'], barcode, script_dir)

        for target_id, gene_name in params["name_dt"].items():
            
            # Load VCF (probably need a method argument)
            vcf_path = f"{params['barcodes_dir']}/{barcode}/call-{method}/reads.target.{gene_name}.vcf"
            vcf_df = load_vcf_as_df(vcf_path)
            
            # Annotate mutations
            vcf_df.insert(0, "mutation", [
                gene.query_for_mutation(
                    chrom=row["chrom"], pos=row["pos"], ref=row["ref"], alt=row["alt"])
                for _, row in vcf_df.iterrows()
            ])

            # Annotate with barcode and gene
            vcf_df.insert(0, "barcode", barcode)
            vcf_df.insert(1, "target_id", target_id)
            vcf_df.insert(2, "gene_name", gene_name)

            # Intersect with mutation frame
            # - Exactly what I return here really depends on what I want to do with the resultant information
            # - E.g., do I want to compute FDR / FNR
            # - Do I want to make plots over all mutations?

            # Store
            dfs.append(vcf_df)


    
    print_footer(t0)

    


