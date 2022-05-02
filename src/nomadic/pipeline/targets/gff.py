import pandas as pd
from dataclasses import dataclass


# ================================================================
# Load and handle gff files
# TODO: this *maybe* should be moved into `lib`
# ================================================================


def load_gff(gff_path):
    """
    Load a Gene Feature Format file from a `gff_path` into
    a pandas data frame

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
        phase: int
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

    return pd.DataFrame(entries)


def extract_gff_attribute(gff_df, extract):
    """Extract attributes from a gff loaded as a data frame"""

    results = []
    for _, row in gff_df.iterrows():
        fields = [r.split("=") for r in row["attributes"].split(";")]
        result = [val for key, val in fields if key == extract]
        results.append(result[0] if result else None)

    return results
