# Quick script to merge BAMs across sequencing experiments
# that are part of the same study
# JHendry, 2024/05/24


import os
import sys
import shutil
import subprocess
import pandas as pd
from dataclasses import dataclass
from nomadic.lib.dirs import ExperimentDirectories
from nomadic.lib.process_bams import samtools_index, samtools_merge


# ================================================================================ #
# Settings
#
# ================================================================================ #


expt_names = {
    "BATCH1": "2024-03-07_ARSUNA-BATCH1-EMMA",
    "BATCH2": "2024-03-13_ARSUNA-BATCH2-EMMA",
    "BATCH3": "2024-05-09_ARSUNA-BATCH3-EMMA",
}


# ================================================================================ #
# Functions
#
# ================================================================================ #


def get_barcodes_dir(expt_name: str) -> str:
    return f"experiments/{expt_name}/nomadic/dorado/sup/single_end_strict/barcodes"


def get_barcode_dir(expt_name: str, barcode: str) -> str:
    barcodes_dir = get_barcodes_dir(expt_name)
    return f"{barcodes_dir}/{barcode}"


def get_pf_bam(expt_name: str, barcode: str) -> str:
    bam_path = f"{get_barcode_dir(expt_name, barcode)}/bams/{barcode}.Pf3D7.final.sorted.bam"
    if not os.path.exists(bam_path):
        raise ValueError(f"Bam not found at {bam_path}!")

    return bam_path


@dataclass
class MetadataRow:
    barcode: str
    sample_id: str
    prior_batch: str = None
    prior_barcode: str = None

    origin: str = None
    is_negative: str = None
    replicate: str = None
    parasitemia: float = None
    parasitemia_log10: float = None

    @classmethod
    def from_dataframe(cls, sample_id: str, new_barcode: str, df: pd.Series):

        TRANSFER_DIRECT = [
            "origin",
            "is_negative",
            "replicate",
            "parasitemia",
            "parasitemia_log10",
        ]

        if df.shape[0] == 1:
            info = df.squeeze()
            return cls(
                sample_id=sample_id,
                barcode=new_barcode,
                prior_batch=info["batch"],
                prior_barcode=info["barcode"],
                **info[TRANSFER_DIRECT],
            )

        # single_get = lambda c: df.squeeze()[c]
        info = df.iloc[1]
        multi_get = lambda c: "|".join(df[c].tolist())

        return cls(
            sample_id=sample_id,
            barcode=new_barcode,
            prior_batch=multi_get("batch"),
            prior_barcode=multi_get("barcode"),
            **info[TRANSFER_DIRECT],
        )


def main(new_expt_name: str, input_metadata: str) -> None:

    # Load
    metadata = pd.read_csv(input_metadata)

    # Create new directories
    new_expt = ExperimentDirectories(new_expt_name)

    # Iterate and merge
    results = []
    for ix, (sample_id, sample_df) in enumerate(metadata.groupby("sample_id")):

        print("-" * 80)
        # Create new barcode directory
        new_barcode = f"barcode{ix:04d}"
        new_barcode_dir = new_expt.get_barcode_dir(new_barcode)
        new_bam_path = f"{new_barcode_dir}/bams/{new_barcode}.Pf3D7.final.sorted.bam"

        old_bams = [
            get_pf_bam(expt_names[batch], barcode)
            for batch, barcode in zip(sample_df["batch"], sample_df["barcode"])
        ]
        n_bams = len(old_bams)

        print(f"Sample ID: {sample_id}")
        print(f"New barcode: {new_barcode}")
        print(f"New directory: {new_barcode_dir}")
        print(f"Num. old bams: {n_bams}")

        if n_bams == 1:
            print("  copying...")
            # Copy operation
            shutil.copy(old_bams[0], new_bam_path)
        else:
            print("  merging...")
            # Merge operation
            unsorted_bam = new_bam_path.replace("sorted.bam", ".bam")
            samtools_merge(old_bams, unsorted_bam)
            cmd = f"samtools sort -Ob -o {new_bam_path} {unsorted_bam}"
            subprocess.run(cmd, check=True, shell=True)

        # Index
        print("  indexing...")
        samtools_index(input_bam=new_bam_path)

        print("Updating metadata record...")
        results.append(MetadataRow.from_dataframe(sample_id, new_barcode, sample_df))

        print("Done.")
        print("-" * 80)

    print("Writing new metadata:")
    new_metadata_csv = f"{new_expt.metadata_dir}/sample_info.csv"
    result_df = pd.DataFrame(results)
    result_df.to_csv(new_metadata_csv, index=False)
    print("Done.")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

    