import os
from nomadic.lib.parsing import build_parameter_dict


# CALLING_FILES = [""]
CFHAPPY_FILES = ["all_targets.sorted.gz.vcf.gz", "all_targets.sorted.gz.vcf.gz"]


def main(epxt_dir, configs, barcodes, method):
    """
    Count for expected target files in a downsampling analysis

    """
    params = build_parameter_dict(expt_dir, config)
    for barcode in barcodes:
        print("=" * 80)
        print(f"Checking barcode: {barcode}")
        print("-"*80)
        barcode_dir = f"{params['barcodes_dir']}/{barcode}"

        # Check variant calling results
        calling_dir = f"{barcode_dir}/calling/{method}/downsample"
        calling_files = os.listdir(calling_dir)

        print(f"Calling directory: {calling_dir}")
        print(f"  Total files: {len(calling_files)}")

        # CFHAPPY files
        cfhappy_dir = f"{barcode_dir}/cfhappy/{method}"
        cfhappy_files = os.listdir(cfhappy_dir)
        print(f"cfhappy directory: {cfhappy_dir}")
        print(f"  Total files: {len(cfhappy_files)}")
    print("="*80)


if __name__ == "__main__":

    # Settings
    expt_dir = "experiments/2021-11-14_strain-validation-flongle-lfb"
    config = "configs/guppy-sup-se.ini"
    barcodes = [f"barcode{i:02d}" for i in range(2, 5)]
    method = "clair3sing"

    # Run
    main(expt_dir, config, barcodes, method)
