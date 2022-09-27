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


# Experiment
EXPT_DIR = "experiments/2022-02-16_zmb-discards-8plex"
CONFIG = "configs/guppy-hac-se.ini"
BARCODES = [f"barcode{i:02d}" for i in range(1, 48)]
METHODS = ["bcftools", "clair3sing", "longshot"]
OUTPUT_DIR = "summary_tables"

# Filtering
KEEP_COLUMNS = [
    "barcode", "method", "Type", "Filter", "Subset",
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


def main():
    """
    Merge HAPPY results and create summary CSVs

    """
    # Load parameters
    params = build_parameter_dict(EXPT_DIR, CONFIG)

    # Iterate over barcodes and methods
    dfs = []
    for barcode, method in product(BARCODES, METHODS):
        
        # Define input directory
        happy_dir = f"{params['barcodes_dir']}/{barcode}/cfhappy/{method}"
        if not os.path.exists(happy_dir):
            print(f"No hap.py comparison for {barcode} and {method}. Skipping.")
            continue

        # Create file matching pattern
        #prefix = f"{barcode}.n([0-9]{'{4}'}).r([0-9]{'{3}'})"
        prefix = "reads"
        suffix = ".all_targets.sorted.gz.extended.csv"
        target_pattern = prefix + suffix
        print(f"Target file pattern: {target_pattern}")
        
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
            # df.insert(2, "replicate", int(file_match.group(2)))
            # df.insert(3, "n_reads", int(file_match.group(1)))
            
            # Clean
            df = df[KEEP_COLUMNS]
            df.query(QUERY, inplace=True)
            
            # Store
            dfs.append(df)
    
    # Combine
    merged_df = pd.concat(dfs)
    merged_df.sort_values(
        ["barcode", "method"]#, "n_reads", "replicate", "Subset"]
    )

    # Write
    output_dir = produce_dir(OUTPUT_DIR, os.path.basename(EXPT_DIR))
    merged_df.to_csv(
        f"{output_dir}/complete_merged.snps.csv",
        index=False
    )
    

if __name__ == "__main__":
    main()