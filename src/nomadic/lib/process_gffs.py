import pandas as pd
from dataclasses import dataclass


# ================================================================
# Load and handle gff files
# 
# ================================================================


def load_gff(gff_path):
    """
    Load a Gene Feature Format file from a `gff_path` into
    a pandas data frame

    param:
        gff_path : str
            Path to .gff file.
    returns:
        _ : Pandas DataFrame
        A .gff file loaded as a pandas dataframe.
        
    """

    # Define an entry in the gff
    @dataclass
    class gffEntry:
        seqid: str
        source: str
        feature: str
        start: int
        end: int
        score: str
        strand: str
        phase: str
        attributes: str

    # Iterate over rows in gff, load entries
    with open(gff_path, "r") as gff:
        entries = []
        for record in gff:
            if record.startswith("#"):
                continue
            fields = record.strip().split("\t")
            entry = gffEntry(*fields)
            entries.append(entry)
    
    # Enforce data types
    gff = pd.DataFrame(entries)
    gff["start"] = gff["start"].astype("int")
    gff["end"] = gff["end"].astype("int")

    return gff


def extract_gff_attribute(gff_df, extract):
    """Extract attributes from a gff loaded as a data frame"""

    results = []
    for _, row in gff_df.iterrows():
        fields = [r.split("=") for r in row["attributes"].split(";")]
        result = [val for key, val in fields if key == extract]
        results.append(result[0] if result else None)

    return results


def add_gff_fields(gff_df, fields=["ID", "Name", "Parent"]):
    
    for field in fields:
        gff_df[field] = extract_gff_attribute(gff_df, extract=field)
    
    return gff_df
