import os
import json
import pandas as pd
from nomadic.lib.parsing import build_parameter_dict  

# Parameters
METADATA_DIR = "experiments/2022-02-16_zmb-discards-8plex/metadata"
MAPPING_CSV = f"{METADATA_DIR}/sample_info_cf_illumina.csv"
OUTPUT_JSON = f"{os.path.dirname(MAPPING_CSV)}/truth_mapping.json"
ILLUMINA_EXPT = "experiments/2022-03-16_ogc-illumina"
ILLUMINA_CONFIG = "configs/default-illumina.ini"
ILLUMINA_CALL_METHOD = "gatk"
TARGET_FILE = "reads.all_targets.sorted.vcf.gz"  # output by ILLUMINA_CALL_METHOD

def main():
    """
    Create a JSON file mapping illumina and
    nanopore barcodes, for clinical discard samples
    
    """

    print("PARAMETERS")
    print(f"  CSV defining mapping between barcoes: {MAPPING_CSV}")
    print(f"  Illumina experiment: {ILLUMINA_EXPT}")
    print(f"  Calling method as truth: {ILLUMINA_CALL_METHOD}")
    print(f"  Target file fomart expected: {TARGET_FILE}")
    print("Done.\n")

    ill_params = build_parameter_dict(
        expt_dir=ILLUMINA_EXPT,
        config=ILLUMINA_CONFIG
    )

    print(f"Loading barcode mapping dataframe...")
    mapping_df = pd.read_csv(MAPPING_CSV)
    print("Done.\n")

    # How are NA handled?
    barcode_to_barcode_map = dict(
        zip(
            mapping_df["barcode"],
            mapping_df["barcode_illumina"].astype("str") # *
        )
    )
    # *Any column with missing data gets coerced to float by pandas

    print("Iterating over nanopore barcodes...")
    mapping_dt = {}
    for ont_barcode, ill_barcode in barcode_to_barcode_map.items():
        
        if ill_barcode == "nan": # Quite ugly, but so it goes with pandas
            # TODO: include controls
            print(f"  No corresponding barcode for {ont_barcode}.")
            continue

        # Get calling directory for Illumina
        ill_barcode_dir = f"{ill_params['barcodes_dir']}/{ill_barcode}"
        ill_call_dir = f"{ill_barcode_dir}/calling/{ILLUMINA_CALL_METHOD}"

        # Check for appropriate file
        found_file = [
            f"{ill_call_dir}/{f}"
            for f in os.listdir(ill_call_dir)
            if f == TARGET_FILE
        ]

        # Add to dictionary
        if not found_file:
            print(f"  No file found at: {ill_call_dir}. Skipping.")
            continue

        mapping_dt[ont_barcode] = found_file[0]
    print("Done.\n")
    
    json.dump(mapping_dt, open(OUTPUT_JSON, "w"))
    print(f"Found {len(mapping_dt)} matching VCFs.")
    print(f"Written to: {OUTPUT_JSON}")
    print("Done.\n")


if __name__ == "__main__":
    main()


        












