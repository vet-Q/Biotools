import subprocess

from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.generic import produce_dir, print_header, print_footer
from ..trim.targets import TARGET_COLLECTION


# --------------------------------------------------------------------------------
# (4) Create pairwise PAF using minimap2
# - Would be good to implement strategy pattern here
# --------------------------------------------------------------------------------


def create_overlap_paf(input_fastq, output_paf):
    """
    Create a PAF of overlaps between reads inside an
    `input_fastq`
    
    """

    # Compute overlap
    cmd = "minimap2 -x ava-ont"
    cmd += f" {input_fastq}"
    cmd += f" {input_fastq}"
    cmd += f" > {output_paf}"

    subprocess.run(cmd, check=True, shell=True)


# --------------------------------------------------------------------------------
# Main script
# 
# --------------------------------------------------------------------------------


def main(expt_dir, config, barcode, target_gene):
    """
    Use `minimap2` to look for overlaps between a set
    of reads deriving from a single `fastq`
    
    """

    # PARSE INPUTS
    script_descrip = "NOMADIC: Look for overlaps between reads from a target gene"
    t0 = print_header(script_descrip)
    script_dir = "coi"
    params = build_parameter_dict(expt_dir, config, barcode)

    target = TARGET_COLLECTION[target_gene]
    print("User inputs:")
    print(f"  Target: {target.name}")
    print(f"  Chrom: {target.chrom}")
    print(f"  Start: {target.start}")
    print(f"  End: {target.end}")
    print("Done.\n")

    # Focus on a single barcode, if specified
    if "focus_barcode" in params:
        params["barcodes"] = [params["focus_barcode"]]

    # ITERATE
    print("Iterating over barcodes...")
    for barcode in params["barcodes"]:
        print("." * 80)
        print(f"Barcode: {barcode}")
        print("." * 80)

        # DIRECTORIES
        coi_dir = produce_dir(params["barcodes_dir"], barcode, script_dir)
        fastq_dir = f"{coi_dir}/fastq_clipped"
        fastq_path = f"{fastq_dir}/reads.target.{target_gene}.clipped.fastq"

        overlap_dir = produce_dir(coi_dir, "overlap")
        paf_path = f"{overlap_dir}/reads.overlap.{target_gene}.paf"
        create_overlap_paf(
            input_fastq=fastq_path,
            output_paf=paf_path
        )

        print("Done.\n")

    print_footer(t0)