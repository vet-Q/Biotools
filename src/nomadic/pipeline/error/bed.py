from abc import ABC, abstractmethod
from dataclasses import dataclass

# TODO:
# - Complete other methods
# - Add a dictionary holding builders


@dataclass
class BEDRow:
    """Encapsulate a single row in a BED file"""

    seqid: str
    start: int
    end: int
    metadata: str = ""

    def as_string(self):
        return f"{self.seqid}\t{self.start}\t{self.end}\t{self.metadata}\n"


class TargetBEDBuilder(ABC):
    def __init__(self, gff_df):
        self.gff_df = gff_df
        self.bed_rows = None

    def set_target(self, target_id):
        self.target_id = target_id

    @abstractmethod
    def find_bed_rows(self):
        """Get information about rows in a bed file"""
        pass

    def write_bed(self, bed_path):
        """
        Write a bed file from information in `self.bed_rows` to
        a file at `bed_path`

        """
        if self.bed_rows is None:
            raise ValueError("No BED rows found. Run `.find_bed_rows()` method first.")

        with open(bed_path, "w") as bed:
            for r in self.bed_rows:
                bed.write(r.as_string())


class ByProteinCodingGene(TargetBEDBuilder):
    """
    Create BED based on `protein_coding_gene` field

    """

    def find_bed_rows(self):
        # Reduce to only target ID
        result_df = (
            self.gff_df.query("feature == 'protein_coding_gene'")
            .query(f"ID == '{self.target_id}'")
            .reset_index()
        )

        # Ensure we have found a single target
        if result_df.shape[0] != 1:
            raise ValueError(f"Can't find information about {self.target_id}")

        # Convert to dictionary
        result_dt = result_df.iloc[0].to_dict()

        # Store
        self.bed_rows = [
            BEDRow(
                seqid=result_dt["seqid"],
                start=int(result_dt["start"]),
                end=int(result_dt["end"]),
                metadata=result_dt["Name"],
            )
        ]


class BetweenStartAndStop(TargetBEDBuilder):
    """
    Between start ad stop

    """
    pass


class AllCDS(TargetBEDBuilder):
    pass
