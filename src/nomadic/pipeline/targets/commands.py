import os
import pandas as pd

import click
from nomadic.pipeline.cli import experiment_options, barcode_option

from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.references import PlasmodiumFalciparum3D7
from nomadic.lib.process_bams import (
    samtools_view,
    samtools_index,
    bedtools_intersect,
    summarise_bam_stats,
)
from .extraction import TargetFactory, write_bed_from_targets
from .plot import plot_dataframe_heat, BalancePlotter
from nomadic.pipeline.qcbams.commands import combine_barcode_dataframes


@click.command(short_help="Analyse amplicon targets.")
@experiment_options
@barcode_option
@click.option(
    "--overview", is_flag=True, help="Produce an overview across all barcodes."
)
def targets(expt_dir, config, barcode, overview):
    """
    Analyse reads overlapping a specific set of targets

    """
    if overview:
        target_extraction_overview(expt_dir, config)
    else:
        target_extraction(expt_dir, config, barcode)


def target_extraction_overview(expt_dir, config):
    """
    Produce plots on the outputs of `target_extraction()`
    across all of the barcodes in an experiment

    """
    # PARSE INPUTS
    script_descrip = "NOMADIC: Overview of target extraction."
    t0 = print_header(script_descrip)
    script_dir = "target-extraction"
    params = build_parameter_dict(expt_dir, config, barcode=None)
    output_dir = produce_dir(params["nomadic_dir"], script_dir)

    # LOAD DATA
    joint_df = combine_barcode_dataframes(
        csv_filename="table.extraction.summary.csv",
        script_dir="target-extraction",
        params=params,
    )

    # Merge with metadata
    merged_df = pd.merge(left=joint_df, right=params["metadata"], on="barcode")

    # Write
    merged_df.to_csv(f"{output_dir}/table.target_coverage.overview.csv", index=False)

    # Iterate over overlap types and plot
    overlap_types = ["complete", "any"]
    for overlap_type in overlap_types:

        # Reduce to plotting information
        plot_df = merged_df.query("overlap == @overlap_type")

        # Pivot for heatmap
        pivot_df = plot_df.pivot(
            index="sample_id", columns="gene_name", values="reads_mapped"
        )

        # Plot heatmap
        plot_dataframe_heat(
            pivot_df,
            cmap="Wistia",
            output_path=f"{output_dir}/plot.reads_mapped.heatmap.{overlap_type}.pdf",
        )

        # Plot stripplot
        plotter = BalancePlotter(
            sample_ids=plot_df["sample_id"],
            gene_names=plot_df["gene_name"],
            values=plot_df["reads_mapped"],
        )
        plotter.set_gene_color_pal("Spectral")
        plotter.plot(
            output_path=f"{output_dir}/plot.reads_mapped.stripplot.{overlap_type}.pdf"
        )

    print_footer(t0)


def target_extraction(expt_dir, config, barcode):
    """
    Extract reads overlapping a set of target genes,
    store summary statistics and write as new .bam files

    """

    # PARSE INPUTS
    script_descrip = "NOMADIC: Extract reads overlapping target regions"
    script_dir = "target-extraction"
    t0 = print_header(script_descrip)
    params = build_parameter_dict(expt_dir, config, barcode)

    # Focus on a single barcode, if specified
    if "focus_barcode" in params:
        params["barcodes"] = [params["focus_barcode"]]

    # Extract reference
    reference = PlasmodiumFalciparum3D7()

    # PREPARE TARGETS
    target_factory = TargetFactory.from_gff_path(reference.gff_path)
    targets = target_factory.get_targets(target_ids=params["target_ids"])
    for target in targets:
        target.name = params["name_dt"][target.ID]

    # ITERATE
    print("Iterating over barcodes...")
    for barcode in params["barcodes"]:
        print("." * 80)
        print(f"Barcode: {barcode}")
        print("." * 80)

        # Prepare to compute per-barcode results
        barcode_results = []

        # Create output directory
        barcode_dir = f"{params['barcodes_dir']}/{barcode}"
        input_dir = f"{barcode_dir}/bams"
        output_dir = produce_dir(barcode_dir, script_dir)

        # Define input bam
        input_bam_path = f"{input_dir}/{barcode}.{reference.name}.final.sorted.bam"

        # Get mapped reads
        mapped_bam_path = f"{output_dir}/reads.mapped.bam"
        samtools_view(input_bam_path, "-F 0x904", mapped_bam_path)
        samtools_index(mapped_bam_path)

        print("  Iterating over targets...")
        for target in targets:
            print(f"\t{target.ID}\t{target.name}")

            target_bam_path = f"{output_dir}/reads.target.{target.name}.bam"

            print(f"\tFinding any overlapping reads...")
            # Any overlap >=2000bp --> removed this to allow for short read overlapping
            samtools_view(mapped_bam_path, f"{target.region}", target_bam_path) # "-m 2000"; can filter by size here
            samtools_index(target_bam_path)

            # Compute summary statistics
            dt = summarise_bam_stats(target_bam_path)
            dt.update(
                {
                    "barcode": barcode,
                    "gene_id": target.ID,
                    "gene_name": target.name,
                    "overlap": "any",
                }
            )
            barcode_results.append(dt)

            # Write a temporary .bed file
            temp_bed_path = target_bam_path.replace(".bam", ".bed")
            write_bed_from_targets(targets=[target], bed_path=temp_bed_path)

            # Find completely overlapping reads
            print(f"\tFinding completely overlapping reads...")
            complete_bam_path = target_bam_path.replace(".bam", ".complete.bam")
            bedtools_intersect(
                input_a=target_bam_path,
                input_b=temp_bed_path,
                args="-F 1.0",
                output=complete_bam_path,
            )
            samtools_index(complete_bam_path)

            # Remove .bed
            os.remove(temp_bed_path)

            # Compute summary statistics
            dt = summarise_bam_stats(complete_bam_path)
            dt.update(
                {
                    "barcode": barcode,
                    "gene_id": target.ID,
                    "gene_name": target.name,
                    "overlap": "complete",
                }
            )
            barcode_results.append(dt)
            print("    Done.")
            print("")

        # Write barcode summary
        print("Writing summary table...")
        barcode_output_path = f"{output_dir}/table.extraction.summary.csv"
        print(f"  to: {barcode_output_path}")
        pd.DataFrame(barcode_results).to_csv(barcode_output_path, index=False)
        print("  Done.")
        print("")
    print("Done.")
    print("")
    print_footer(t0)
