import os
import re
import pandas as pd
from itertools import product
from nomadic.lib.generic import produce_dir
from nomadic.lib.parsing import build_parameter_dict


# ================================================================================
# Parameters
#
# ================================================================================


# Output
OUTPUT_DIR = "summary_tables"

# Filtering
KEEP_COLUMNS = [
    "barcode", "method", "n_reads", "replicate", 
    "Type", "Filter", "Subset",
    "METRIC.Recall", "METRIC.Precision", "METRIC.F1_Score", "METRIC.Frac_NA",
    "TRUTH.TOTAL", "TRUTH.TP", "TRUTH.FN",
    "QUERY.TOTAL", "QUERY.TP", "QUERY.FP", "QUERY.UNK",
    "FP.gt", "FP.al",
    "Subset.IS_CONF.Size"
]
QUERY = "Filter == 'PASS' and Type == 'SNP'"


# ================================================================================
# Script
#
# ================================================================================


def main(expt_dir, config, barcodes, methods):
    """
    Merge HAPPY results and create summary CSVs

    """
    # Load parameters
    params = build_parameter_dict(expt_dir, config)
    expt_name = os.path.basename(expt_dir)

    # Iterate over barcodes and methods
    dfs = []
    for barcode, method in product(barcodes, methods):
        print("Aggregating:")
        print(f"  Barcode {barcode}")
        print(f"  Method: {method}")
        
        # Define input directory
        happy_dir = f"{params['barcodes_dir']}/{barcode}/cfhappy/{method}"
        if not os.path.exists(happy_dir):
            print(f"No hap.py comparison for {barcode} and {method}. Skipping.")
            continue

        # Create file matching pattern
        prefix = f"{barcode}.n([0-9]{'{4}'}).r([0-9]{'{3}'})"
        #prefix = "reads"
        suffix = ".all_targets.sorted.gz.extended.csv"
        target_pattern = prefix + suffix
        #print(f"Target file pattern: {target_pattern}")
        
        for file in os.listdir(happy_dir):
            
            # Check for match
            file_match = re.match(target_pattern, file)
            
            # Skip if not
            if not file_match:
                continue
                
            # Load
            df = pd.read_csv(f"{happy_dir}/{file}")
            
            # Annotate
            df.insert(0, "barcode", barcode)
            df.insert(1, "method", method)
            df.insert(2, "replicate", int(file_match.group(2)))
            df.insert(3, "n_reads", int(file_match.group(1)))
            
            # Clean
            df = df[KEEP_COLUMNS]
            df.query(QUERY, inplace=True)
            
            # Store
            dfs.append(df)
        print(f"  Number of aggregated samples: {len(dfs)}")
        print("Done.")
    
    merged_df = pd.concat(dfs)
    merged_df.sort_values(
        ["barcode", "method"]#, "n_reads", "replicate", "Subset"]
    )

    # Write
    print("Writing output...")
    output_dir = produce_dir(OUTPUT_DIR, os.path.basename(EXPT_DIR))
    output_path = f"{output_dir}/sup/fig4_aggregated-summary.csv"
    print(f"  to: {output_path}")
    merged_df.to_csv(
        output_path,
        index=False
    )
    print("Done.")
    

if __name__ == "__main__":

    # Selecting
    expt = "flongle-lfb-sup"

    # Clinical discards
    if expt == "clinical":
        EXPT_DIR = "experiments/2022-02-16_zmb-discards-8plex"
        CONFIG = "configs/guppy-hac-se.ini"
        BARCODES = [f"barcode{i:02d}" for i in range(1, 48)]
        METHODS = ["bcftools", "clair3sing", "longshot"]
    elif expt == "flongle-lfb":
        EXPT_DIR = "experiments/2021-11-14_strain-validation-flongle-lfb"
        CONFIG = "configs/guppy-hac-se.ini"
        BARCODES = [f"barcode{i:02d}" for i in range(2, 5)]
        #METHODS = ["bcftools", "clair3sing", "longshot"]
        METHODS = ["clair3sing"]
    elif expt == "flongle-lfb-sup":
        EXPT_DIR = "experiments/2021-11-14_strain-validation-flongle-lfb"
        CONFIG = "configs/guppy-sup-se.ini"
        BARCODES = [f"barcode{i:02d}" for i in range(2, 5)]
        #METHODS = ["bcftools", "clair3sing", "longshot"]
        METHODS = ["clair3sing"]

    # Run
    main(
        expt_dir=EXPT_DIR,
        config=CONFIG,
        barcodes=BARCODES,
        methods=METHODS
    )