import os
import warnings
from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.references import PlasmodiumFalciparum3D7
from .callers import caller_collection
from .annotator import VariantAnnotator
from .merger import VariantMerger


# from nomadic.lib.process_vcfs import bcftools_view, bcftools_concat, bcftools_sort, bcftools_index

def quickcall(expt_dir: str, config: str, barcode: str, bed_path: str, method: str, overview: bool=False) -> None:
    if overview:
        quickcall_merge(expt_dir, config, bed_path, method)
    else:
        quickcall_single(expt_dir, config, barcode, bed_path, method)



def quickcall_single(expt_dir: str, config: str, barcode: str, bed_path: str, method: str) -> None:
    """
    Effectively an improved approach to variant calling
    vs. the old `call`

    NB:
    - Now we are *not* necessarily filtering the BAM file
    of chimeric or secondary reads
    - This may occur by the variant caller, but not ensured
    - Should consider introducing a filtering step for the bam file,
    probably via `samtools view`

    """

    # PARSE INPUTS
    script_descrip = "NOMADIC: Call variants without downsampling"
    t0 = print_header(script_descrip)
    script_dir = f"quickcall/{method}"  # this will be easier to iterate overs
    params = build_parameter_dict(expt_dir, config, barcode)

    # Define reference genome
    reference = PlasmodiumFalciparum3D7()

    # Focus on a single barcode, if specified
    if "focus_barcode" in params:
        params["barcodes"] = [params["focus_barcode"]]

    # Iterate over barcodes
    for barcode in params["barcodes"]:
        # Define input and output directory
        print("-" * 80)
        print(f"Running {method} for: {barcode}")
        barcode_dir = f"{params['barcodes_dir']}/{barcode}"
        input_dir = f"{barcode_dir}/bams"
        output_dir = produce_dir(barcode_dir, script_dir)

        # Path to *complete* bam file
        bam_path = f"{input_dir}/{barcode}.{reference.name}.final.sorted.bam"
        vcf_path = f"{output_dir}/{barcode}.{reference.name}.{method}.unfiltered.vcf.gz"

        # TODO:
        # -> I *could* and probably *should* filter at this point

        # Get variant caller and call
        CallingMethod = caller_collection[method]
        caller = CallingMethod(fasta_path=reference.fasta_path)

        print("Calling variants...")
        caller.run(bam_path, vcf_path, sample_name=barcode)
        print("Done.")
        print("")

        print("Filtering variants...")
        filtered_vcf = vcf_path.replace(".unfiltered.vcf.gz", ".filtered.vcf.gz")
        caller.filter(output_vcf=filtered_vcf, bed_path=bed_path)

        print("Reducing to biallelic...")
        biallelic_vcf = vcf_path.replace(".unfiltered.vcf.gz", ".biallelic.filtered.vcf.gz")
        caller.filter(output_vcf=biallelic_vcf, to_biallelic=True, bed_path=bed_path)

        # Annotation
        print("Annotating variants...")
        annotator = VariantAnnotator(
            biallelic_vcf, # Note that we annotated only filtered VCF
            bed_path,
            reference,
            output_dir=output_dir
        )
        annotator.run()
        annotator.convert_to_tsv()
        print("Done.")
        print("")

    print_footer(t0)



def quickcall_merge(expt_dir: str, config: str, bed_path: str, method: str) -> None:
    """
    Merge the VCF files and write a final TSV
    
    """
    # PARSE INPUTS
    script_descrip = "NOMADIC: Merge quick variant calls across experiment"
    t0 = print_header(script_descrip)
    params = build_parameter_dict(expt_dir, config)
    script_dir = f"quickcall/{method}"  # this will be easier to iterate overs
    output_dir = produce_dir(params["nomadic_dir"], script_dir)

    # Define reference genome
    reference = PlasmodiumFalciparum3D7()

    # Iterate over barcodes
    vcfs = []
    for barcode in params["barcodes"]:
        # Define input and output directory
        print("-" * 80)
        print(f"Running {method} for: {barcode}")
        barcode_dir = f"{params['barcodes_dir']}/{barcode}"
        input_dir = f"{barcode_dir}/quickcall/{method}"

        # Path to *complete* bam file
        vcf_path = f"{input_dir}/{barcode}.{reference.name}.{method}.filtered.vcf.gz"

        if os.path.exists(vcf_path):
            vcfs.append(vcf_path)
        else:
            warnings.warn(f"No VCF file found at {vcf_path}! Skipping.")

    print(f"Found {len(vcfs)} total.")

    # Merging
    output_vcf = f"{output_dir}/merged.vcf.gz"
    merger = VariantMerger(vcfs)
    merger.run(output_vcf)

    # Annotation
    print("Annotating variants...")
    annotator = VariantAnnotator(
        output_vcf, # Note that we annotated only filtered VCF
        bed_path,
        reference,
        output_dir=output_dir
    )
    annotator.run()
    annotator.convert_to_tsv()
    print("Done.")
    print("")

    print_footer(t0)
