import os
import pandas as pd
import numpy as np
from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.references import PlasmodiumFalciparum3D7
from nomadic.lib.process_bams import samtools_depth, bedtools_intersect
from nomadic.lib.process_gffs import write_gff_to_bed


def depth(expt_dir, config, barcode, amplicon_bed_path, cds_bed_path):
    """
    Analyse depth profiles across a set of amplicons

    """
    # Parse inputs
    script_descrip = "NOMADIC: Compute read depth over amplicon regions"
    script_dir = "depth"
    t0 = print_header(script_descrip)
    print("Parse inputs...")
    print(f"  Experiment dir.: {expt_dir}")
    print(f"  Configuration path: {config}")
    print(f"  Amplicon BED file: {amplicon_bed_path}")
    print(f"  CDS BED file: {cds_bed_path}")
    print("Done.")
    print("")

    # Load parameters
    params = build_parameter_dict(expt_dir=expt_dir, config=config, barcode=barcode)
    reference = PlasmodiumFalciparum3D7()

    # Focus on a single barcode, if specified
    if "focus_barcode" in params:
        params["barcodes"] = [params["focus_barcode"]]

    # Iterate over barcodes
    for barcode in params["barcodes"]:

        print("." * 80)
        print(f"Barcode: {barcode}")
        print("." * 80)

        # Define relevant paths
        barcode_dir = f"{params['barcodes_dir']}/{barcode}"
        bam_file = f"{barcode}.{reference.name}.final.sorted.bam"
        bam_path = f"{barcode_dir}/bams/{bam_file}"
        depth_dir = produce_dir(barcode_dir, script_dir)
        depth_path = f"{depth_dir}/{bam_file.replace('.bam','.depth')}"

        # Compute depths
        print("Computing depths...")
        samtools_depth(
            input_bam=bam_path, region_bed=amplicon_bed_path, output_path=depth_path
        )

        # Load and convert to bed, for merging
        depth_df = pd.read_csv(depth_path, sep="\t", names=["chrom", "start", "depth"])
        depth_df.insert(2, "end", depth_df["start"] + 1)
        depth_bed_path = f"{depth_path}.bed"
        write_gff_to_bed(
            df=depth_df,
            bed_path=depth_bed_path,
            bed_columns=["chrom", "start", "end", "depth"],
        )
        print(f"  Total bases covered: {depth_df.shape[0]}bp")
        print(f"  Average depth: {depth_df['depth'].mean():.2f}")

        # Merge with amplicons
        print("Annotating amplicons...")
        depth_annotate_path = depth_bed_path.replace(".bed", ".amplicon.bed")
        bedtools_intersect(
            input_a=depth_bed_path,
            input_b=amplicon_bed_path,
            args="-wa -wb",
            output=depth_annotate_path,
        )

        # Merge with CDS
        print("Annotating CDS...")
        depth_annotate_cds_path = depth_bed_path.replace(".bed", ".amplicon.cds.bed")
        bedtools_intersect(
            input_a=depth_annotate_path,
            input_b=cds_bed_path,
            args="-wa -wb -loj",
            output=depth_annotate_cds_path,
        )
        # Clean for space
        os.remove(depth_annotate_path)
        os.remove(depth_bed_path)

        # Create summary table
        complete_df = pd.read_csv(
            depth_annotate_cds_path,
            sep="\t",
            names=[
                "chrom",
                "start",
                "end",
                "depth",
                "amp_chrom",
                "amp_start",
                "amp_end",
                "amp_target_id",
                "cds_chrom",
                "cds_start",
                "cds_end",
                "cds_target_id",
            ],
        )
        summary_df = (
            complete_df.groupby("cds_target_id")
            .agg(
                n_bases=pd.NamedAgg("depth", len),
                mean_depth=pd.NamedAgg("depth", np.mean),
                median_depth=pd.NamedAgg("depth", np.median),
                frac_below_1X=pd.NamedAgg("depth", lambda x: sum(x <= 0) / len(x)),
                frac_below_10X=pd.NamedAgg("depth", lambda x: sum(x <= 10) / len(x)),
                frac_below_20X=pd.NamedAgg("depth", lambda x: sum(x <= 20) / len(x)),
            )
            .reset_index()
        )
        summary_df.insert(0, "barcode", barcode)
        summary_path = f"{depth_dir}/table.{barcode}.depth_summary.csv"
        summary_df.to_csv(summary_path, index=False)

        print("Outputs written to:")
        print(f" {depth_path}")
        print(f" {depth_annotate_cds_path}")
        print(f"  {summary_path}")
        print("Done.\n")
    print_footer(t0)
