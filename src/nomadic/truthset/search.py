import json
import click
import pandas as pd
pd.options.mode.chained_assignment = None 
from nomadic.lib.generic import produce_dir, print_header, print_footer
from nomadic.lib.references import (
    PlasmodiumFalciparum3D7,
    PlasmodiumFalciparumDd2,
    PlasmodiumFalciparumGB4,
    PlasmodiumFalciparumHB3
)
from nomadic.lib.process_gffs import load_gff, add_gff_fields


# ================================================================
# Parameters
#
# ================================================================


DEFAULT_MAP_PATH = "src/nomadic/truthset/id_to_names.json"


# ================================================================
# Search for genes of interest across trains based on their name
#
# ================================================================


def search_for_target_genes_by_name(gff_df, id_to_name_mapping, add_standard_id=True, verbose=True):
    """
    Reduce a GFF dataframe `gff_df` to rows that correspond to
    gene names defined in an set `id_to_name_mapping`

    ISSUE:
    - If genes missing == genes duplicated; artificially
    looks like you found a result for each gene
    
    param:
        gff_df: DataFrame
            Gene format file loaded as data frame.
        id_to_name_mapping: dict
            Keys are IDs in a reference genome, values are
            possible gene names for those IDs.
            
    returns
        target_df: DataFrame
            Rows corresponding to genes found in `id_to_name_mapping`
    
    """
    
    # Search for individual targets
    target_rows = []
    target_counts = {}
    for target_id, target_names in id_to_name_mapping.items():
        
        r = gff_df.query("name in @target_names")
        
        if r.shape[0]:
            if add_standard_id:
                r["standard_id"] = target_id
            
            # If we append here, means at least 1 target was found
            target_rows.append(r)
            target_counts[target_names[0]] = r.shape[0]
        
        elif verbose:
            print(f" For ID={target_id}, no genes in {', '.join(target_names)} found.")
            target_counts[target_names[0]] = 0
            
    # Combine across targets
    target_df = pd.concat(target_rows)
            
    if verbose:
        print(f"Total of {len(target_rows)}/{len(id_to_name_mapping)} targets found.")
        sep = "\n"
        print(f"{sep.join([f'  {k:>10}: {v:<3}' for k,v in target_counts.items()])}")
        
    return target_df



# ================================================================
# Main script
#
# ================================================================



@click.command(short_help="Search for genes.")
@click.option(
    "-n",
    "--name_path",
    type=str,
    help="Path to .json file mapping IDs to names.", 
    default=DEFAULT_MAP_PATH
)
def search(name_path):
    """
    Search for genes across various strains on the basis of their name
    
    """
    # Create output directory
    t0 = print_header("TRUTHSET: Searching for genes across various P.f. strains, by their name")
    output_dir = produce_dir("resources", "truthsets")

    # Define reference genomes
    references = [
        PlasmodiumFalciparum3D7(),
        PlasmodiumFalciparumDd2(),
        PlasmodiumFalciparumGB4(),
        PlasmodiumFalciparumHB3()
    ]

    # Load mapping
    id_to_names = json.load(open(name_path, "r"))

    # CREATE .GFF ACROSS ALL STRAINS
    # Iterate over references
    target_dfs = []
    for reference in references:
        
        print("." * 80)
        print(f"Reference: {reference.name}")
        print("." * 80)
        
        # PROCESS GFF
        print(f"Processing .gff: {reference.gff_path}")
        gff_df = (
            load_gff(reference.gff_path)
            .query("feature == 'protein_coding_gene'")
            .reset_index(drop=True)
        )
        gff_df = add_gff_fields(gff_df, fields=["ID", "Name"])
        gff_df.rename({"Name": "name"}, axis=1, inplace=True)
        
        # EXTRACT TARGET GENES
        print("Finding target genes...")
        target_df = search_for_target_genes_by_name(gff_df, id_to_names)

        print("Storing...")
        target_df.insert(0, "reference", reference.name)
        target_dfs.append(target_df)
        print("Done.")
        print("")
    print("")

    # Combine and write
    print("Writing...")
    output_path = f"{output_dir}/table.target_genes.csv"
    print(f"  Output file written to: {output_path}")
    combined_df = pd.concat(target_dfs)
    combined_df.to_csv(output_path, index=False)
    print("Done.")
    print("")
    print_footer(t0)





