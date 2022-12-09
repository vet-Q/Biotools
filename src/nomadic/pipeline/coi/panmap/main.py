


from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.generic import produce_dir, print_header, print_footer
from nomadic.lib.references import (
    PlasmodiumFalciparum3D7,
    PlasmodiumFalciparumDd2,
    PlasmodiumFalciparumGB4,
    PlasmodiumFalciparumHB3
)

from nomadic.pipeline.coi.trim.targets import TARGET_COLLECTION
from .mappers import Minimap2PAF


def main(expt_dir, config, barcode, target_gene):
    """
    Filter and trim all reads in a BAM file to overlap a `target_gene` 
    and span `start` and `end` positions; then convert to FASTQ
    
    """

    # PARSE INPUTS
    script_descrip = "NOMADIC: Trim target BAM file for COI analysis"
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

    # DEFINE REFERNCE GENOMES
    references = [
        PlasmodiumFalciparum3D7(),
        PlasmodiumFalciparumDd2(),
        PlasmodiumFalciparumGB4(),
        PlasmodiumFalciparumHB3()
    ]
    print(f"Mapping to {len(references)} references: {', '.join([r.name for r in references])}")

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

        for reference in references:
            print(f"SPECIES: {reference.name}")

            # Produce required directory
            output_dir = produce_dir(coi_dir, "panmap")
            output_bam = f"{output_dir}/{barcode}.{reference.name}.{target_gene}.sorted.paf"

            # Instantiate mapper
            mapper = Minimap2PAF(reference)

            # Map
            print("Mapping...")
            mapper.map_from_fastqs(fastq_dir=fastq_dir)
            mapper.run(output_bam)
            print("Done.\n")
        print("Done.\n")

    print_footer(t0)


