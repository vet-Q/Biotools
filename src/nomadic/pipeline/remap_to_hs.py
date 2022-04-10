import subprocess

from nomadic.lib.generic import print_header, print_footer
from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.process_bams import samtools_index
from nomadic.lib.references import (
    PlasmodiumFalciparum3D7,
    HomoSapiens,
)


# ================================================================
# Run Minimap2, but only on unmapped reads from a .bam
#
# ================================================================


def remap_with_minimap2(input_bam, output_bam, reference, flags="--eqx --MD"):
    """Remap unmapped reads to a new reference genome"""

    cmd = f"samtools view -f 0x004 {input_bam}"  # get unmapped reads
    cmd += " | samtools fastq"  # convert back to .fastq
    cmd += f" | minimap2 -ax map-ont {flags} {reference.fasta_path} -"  # remap
    cmd += " | samtools view -S -b -"  # convert to .bam
    cmd += f" | samtools sort -o {output_bam}"  # sort and write

    subprocess.run(cmd, shell=True, check=True)

    return None


# ================================================================
# Main script, run from `cli.py`
#
# ================================================================


def main(expt_dir, config, barcode):
    """
    Remap all reads that failed to map to P.f. referece genome
    to human referece genome

    """

    # PARSE INPUTS
    script_descrip = "NOMADIC: Re-map unmapped reads to Homo Sapiens"
    t0 = print_header(script_descrip)
    script_dir = "bams"
    params = build_parameter_dict(expt_dir, config, barcode)

    # Focus on a single barcode, if specified
    if "focus_barcode" in params:
        params["barcodes"] = [params["focus_barcode"]]

    # Define reference genomes
    pf_reference = PlasmodiumFalciparum3D7()
    hs_reference = HomoSapiens()

    # ITERATE
    print("Iterating over barcodes...")
    for barcode in params["barcodes"]:
        print("." * 80)
        print(f"Barcode: {barcode}")
        print("." * 80)

        # Define input and output bams
        barcode_dir = f"{params['barcodes_dir']}/{barcode}/{script_dir}"
        input_bam = f"{barcode_dir}/{barcode}.{pf_reference.name}.final.sorted.bam"
        output_bam = f"{barcode_dir}/{barcode}.{hs_reference.name}.final.sorted.bam"

        # Remap
        print("Remapping to H.s...")
        remap_with_minimap2(input_bam, output_bam, hs_reference)
        # Index
        print("Indexing...")
        samtools_index(input_bam=output_bam)
        print("Done.")
        print("")
    print_footer(t0)
