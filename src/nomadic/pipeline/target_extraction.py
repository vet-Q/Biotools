import os
import pandas as pd
from dataclasses import dataclass

from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.references import PlasmodiumFalciparum3D7
from nomadic.lib.process_bams import (
    samtools_view,
    samtools_index,
    bedtools_intersect,
    summarise_bam_stats,
)


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


# ================================================================
# Main script, run from `cli.py`
#
# ================================================================


def main(expt_dir, config, barcode):
    """
    Extract reads overlapping a set of target genes,
    store summary statistics and write as new .bam files

    """

    # PARSE INPUTS
    script_descrip = "NOMADIC: Extract reads overlapping target regions"
    script_dir = "target-extraction"
    t0 = print_header(script_descrip)
    params = build_parameter_dict(expt_dir, config, barcode)

    # Focus on a single barcode, if specified
    if "focus_barcode" in params:
        params["barcodes"] = [params["focus_barcode"]]

    # Extract reference
    reference = PlasmodiumFalciparum3D7()

    # PREPARE TARGETS
    target_factory = TargetFactory.from_gff_path(reference.gff_path)
    targets = target_factory.get_targets(target_ids=params["target_ids"])
    for target in targets:
        target.name = params["name_dt"][target.ID]

    # ITERATE
    print("Iterating over barcodes...")
    for barcode in params["barcodes"]:
        print("." * 80)
        print(f"Barcode: {barcode}")
        print("." * 80)

        # Prepare to compute per-barcode results
        barcode_results = []

        # Create output directory
        barcode_dir = f"{params['barcodes_dir']}/{barcode}"
        input_dir = f"{barcode_dir}/bams"
        output_dir = produce_dir(barcode_dir, script_dir)

        # Define input bam
        input_bam_path = f"{input_dir}/{barcode}.{reference.name}.final.sorted.bam"

        # Get mapped reads
        mapped_bam_path = f"{output_dir}/reads.mapped.bam"
        samtools_view(input_bam_path, "-F 0x904", mapped_bam_path)
        samtools_index(mapped_bam_path)

        print("  Iterating over targets...")
        for target in targets:
            print(f"\t{target.ID}\t{target.name}")

            target_bam_path = f"{output_dir}/reads.target.{target.name}.bam"

            print(f"\tFinding any overlapping reads...")
            # Any overlap >=2000bp
            samtools_view(mapped_bam_path, f"{target.region} -m 2000", target_bam_path)

            # Compute summary statistics
            dt = summarise_bam_stats(target_bam_path)
            dt.update(
                {
                    "barcode": barcode,
                    "gene_id": target.ID,
                    "gene_name": target.name,
                    "overlap": "any",
                }
            )
            barcode_results.append(dt)

            # Write a temporary .bed file
            temp_bed_path = target_bam_path.replace(".bam", ".bed")
            write_bed_from_targets(targets=[target], bed_path=temp_bed_path)

            # Find completely overlapping reads
            print(f"\tFinding completely overlapping reads...")
            complete_bam_path = target_bam_path.replace(".bam", ".complete.bam")
            bedtools_intersect(
                input_a=target_bam_path,
                input_b=temp_bed_path,
                args="-F 1.0",
                output=complete_bam_path,
            )

            # Remove .bed
            os.remove(temp_bed_path)

            # Compute summary statistics
            dt = summarise_bam_stats(complete_bam_path)
            dt.update(
                {
                    "barcode": barcode,
                    "gene_id": target.ID,
                    "gene_name": target.name,
                    "overlap": "complete",
                }
            )
            barcode_results.append(dt)
            print("    Done.")
            print("")

        # Write barcode summary
        print("Writing summary table...")
        barcode_output_path = f"{output_dir}/table.extraction.summary.csv"
        print(f"  to: {barcode_output_path}")
        pd.DataFrame(barcode_results).to_csv(barcode_output_path, index=False)
        print("  Done.")
        print("")
    print("Done.")
    print("")
    print_footer(t0)
