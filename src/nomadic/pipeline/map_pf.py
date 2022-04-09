import os
import sys
from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.process_bams import samtools_index
from nomadic.lib.mapping import Mapper
from nomadic.lib.references import PlasmodiumFalciparum3D7


def main(expt_dir, config, barcode):
    """
    Map .fastq files found in the experiment directory `expt_dir` to the
    P. falciparum reference genome.
    
    """

    # PARSE INPUTS
    script_descrip = "NOMADIC: Map .fastq files to Plasmodium falciparum"
    t0 = print_header(script_descrip)
    script_dir = "bams"
    params = build_parameter_dict(expt_dir, config, barcode)

    # Focus on a single barcode, if specified
    if "focus_barcode" in params:
        params["barcodes"] = [params["focus_barcode"]]

    # Define reference genomes
    references = [
        PlasmodiumFalciparum3D7(),
    ]

    # ITERATE
    print("Iterating over barcodes and references...")
    for barcode in params["barcodes"]:
        print("." * 80)
        print(f"Barcode: {barcode}")
        print("." * 80)

        # Define .fastq path
        fastq_dir = f"{params['fastq_dir']}/{barcode}"
        n_fastqs = len(os.listdir(fastq_dir))
        print(f"Discovered {n_fastqs} .fastq files.")
        if n_fastqs == 0:
            continue

        results = []
        for reference in references:
            print(f"SPECIES: {reference.name}")

            # Produce required directory
            barcode_dir = produce_dir(params["barcodes_dir"], barcode, script_dir)
            output_bam = f"{barcode_dir}/{barcode}.{reference.name}.final.sorted.bam"

            # Instantiate mapper
            mapper = Mapper(fastq_dir, reference)

            # Map
            print("Mapping...")
            mapper.run_minimap2(output_bam)

            # Index
            print("Indexing...")
            samtools_index(input_bam=output_bam)
            print("Done.")
            print("")
        print("")
    print_footer(t0)

    
if __name__ == "__main__":

    # Then I probably just add the decorator? Or...?
    main()