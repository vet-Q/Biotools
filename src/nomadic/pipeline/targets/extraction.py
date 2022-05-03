import pandas as pd
from dataclasses import dataclass
from nomadic.lib.process_gffs import load_gff, extract_gff_attribute


# ================================================================
# Create target region objects
#
# ================================================================


@dataclass
class Target:
    """
    Define a genomic target region

    """

    ID: str
    chrom: str
    start: int
    end: int
    name: str = ""
    length: int = 0
    region: str = ""

    def __post_init__(self):
        self.length = self.end - self.start
        self.region = f"{self.chrom}:{self.start}-{self.end}"

    @classmethod
    def from_gff_record(cls, record):
        """
        A target region can created directly from a `record` within
        a .gff file, if that file has been loaded as a pandas data frame

        """

        return cls(
            ID=record["ID"],
            chrom=record["seqid"],
            start=int(record["start"]),
            end=int(record["end"]),
        )


class TargetFactory:
    def __init__(self, gff_df):
        """
        Create a list of Target objects from a GFF file and
        a set of target IDs

        """
        # Data Frames
        self.gff_df = gff_df
        self.cds_df = self._reduce_to_cds()

    @classmethod
    def from_gff_path(cls, gff_path):
        """
        Initialise from a path to a gff file

        """

        return cls(load_gff(gff_path))

    def _reduce_to_cds(self):
        """
        Restrict the `gff_df` to only coding sequence (CDS)
        features

        """

        cds_df = self.gff_df.query("feature == 'CDS'")
        cds_df.insert(0, "ID", extract_gff_attribute(cds_df, extract="ID"))

        return cds_df

    def get_target(self, target_id):
        """
        Return a Target object for a specific `target_id`

        """

        # Reduce to `target_id` rows
        target_df = self.cds_df.loc[
            [ID.startswith(target_id) for ID in self.cds_df["ID"]]
        ]
        assert target_df.shape[0] > 0, f"No CDS found for {target_id}."

        # Extract chromosome as string
        chrom = target_df["seqid"].unique()[0]

        # Get full extent of cds from across rows
        start = int(target_df["start"].min())
        end = int(target_df["end"].max())

        return Target(ID=target_id, chrom=chrom, start=start, end=end)

    def get_targets(self, target_ids):
        """
        Return a list of Target objects from a list
        of `target_ids`

        """

        targets = []
        for target_id in target_ids:
            targets.append(self.get_target(target_id))

        return targets


def write_bed_from_targets(targets, bed_path):
    """
    Write a .bed file from `targets` a list of Target objects,
    to a provided `bed_path`

    """

    with open(bed_path, "w") as bed:
        for target in targets:
            bed.write(f"{target.chrom}\t{target.start}\t{target.end}\t{target.name}\n")
