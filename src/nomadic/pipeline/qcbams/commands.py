import pandas as pd

pd.options.mode.chained_assignment = None

import click
from nomadic.pipeline.cli import experiment_options, barcode_option
from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.references import (
    PlasmodiumFalciparum3D7,
    HomoSapiens,
)
from .io import load_alignment_information
from .classify import reduce_to_read_dataframe, convert_column_to_ordered_category
from .plot import (
    MappingStatesAndColors,
    HISTOGRAM_STATS,
    JointHistogramPlotter,
    ReadHistogramPlotter,
    barplot_states,
)


# ================================================================
# Load dataframes across all barcodes in the experiment
#
# ================================================================


def combine_barcode_dataframes(csv_filename, script_dir, params):
    """
    Combine a particular dataframe across all barcodes

    """

    dfs = []
    for barcode in params["barcodes"]:
        csv_path = f"{params['barcodes_dir']}/{barcode}/{script_dir}/{csv_filename}"
        df = pd.read_csv(csv_path)
        if not "barcode" in df.columns:
            df.insert(0, "barcode", barcode)
        dfs.append(df)

    return pd.concat(dfs)


# ================================================================
# Main script, run from `cli.py`
#
# ================================================================


@click.command(short_help="QC analysis of .bam files.")
@experiment_options
@barcode_option
@click.option(
    "--overview", is_flag=True, help="Produce an overview across all barcodes."
)
def qcbams(expt_dir, config, barcode, overview):
    """
    Run a quality control analysis of .bam files generated from
    `nomadic map` and `nomadic remap`

    """
    if overview:
        qcbams_overview(expt_dir, config)
    else:
        qcbams_individual(expt_dir, config, barcode)


def qcbams_overview(expt_dir, config):
    """
    Produce a summary of .bam QC across all barcodes within an expeirment

    """
    # PARSE INPUTS
    script_descrip = "NOMADIC: Overview of .bam file QC across all barcodes"
    t0 = print_header(script_descrip)
    script_dir = "qc-bams"
    params = build_parameter_dict(expt_dir, config, barcode=None)
    output_dir = produce_dir(params["nomadic_dir"], script_dir)

    # Instantiate mapping states and colors
    msc = MappingStatesAndColors()

    # Iterate over states and produce barplots
    states = ["primary_state", "secondary_state"]
    for state in states:
        # Load data
        joint_df = combine_barcode_dataframes(
            csv_filename=f"table.size.{state}.csv", script_dir=script_dir, params=params
        )

        # Tall to wide
        joint_df = joint_df.pivot(
            index="barcode", columns="group", values="n_reads"
        ).reset_index()

        # Merge metadata
        merged_df = pd.merge(left=joint_df, right=params["metadata"], on="barcode")
        assert (
            joint_df.shape[0] == merged_df.shape[0]
        ), "Barcodes lost during merge with metadata."

        # Write
        merged_df.to_csv(f"{output_dir}/table.mapping.{state}.csv", index=False)

        # Isolate to key columns for plotting
        plot_df = merged_df[msc.level_sets[state][::-1]]
        plot_df.index = merged_df["sample_id"]

        # Barplot
        barplot_states(
            plot_df,
            colors=msc.color_sets[state][::-1],
            show_per="pf_uniq_mapped" if state == "secondary_state" else "pf_mapped",
            output_path=f"{output_dir}/plot.mapping.{state}.pdf",
        )

        # Now plot histograms
        for histogram_stat in HISTOGRAM_STATS:
            # Load data
            joint_df = combine_barcode_dataframes(
                csv_filename=f"table.bin_counts.{histogram_stat.stat}.{state}.csv",
                script_dir=script_dir,
                params=params,
            )

            # Aggregate
            sum_df = joint_df.groupby(["bin_lower", "bin_higher"]).sum().reset_index()

            # Plot histogram
            plotter = JointHistogramPlotter(
                sum_df, msc.level_sets[state], msc.color_sets[state]
            )
            plotter.set_histogram_stats(**histogram_stat.__dict__)
            plotter.plot_histogram(
                output_path=f"{output_dir}/plot.hist.{histogram_stat.stat}.{state}.pdf"
            )
    print_footer(t0)


def qcbams_individual(expt_dir, config, barcode):
    """
    Create a series of histogram summaries of reads
    within .bam files

    """

    # PARSE INPUTS
    script_descrip = "NOMADIC: Quality control analysis of .bam files"
    t0 = print_header(script_descrip)
    script_dir = "qc-bams"
    params = build_parameter_dict(expt_dir, config, barcode)

    # Focus on a single barcode, if specified
    if "focus_barcode" in params:
        params["barcodes"] = [params["focus_barcode"]]

    # Define reference genomes
    pf_reference = PlasmodiumFalciparum3D7()
    hs_reference = HomoSapiens()

    # Instantiate mapping states and colors
    msc = MappingStatesAndColors()

    # ITERATE
    print("Iterating over barcodes...")
    for barcode in params["barcodes"]:
        print("." * 80)
        print(f"Barcode: {barcode}")
        print("." * 80)

        # Define directories
        barcode_dir = f"{params['barcodes_dir']}/{barcode}"
        bam_dir = f"{barcode_dir}/bams"
        output_dir = produce_dir(barcode_dir, script_dir)

        # Define input bams
        pf_bam_path = f"{bam_dir}/{barcode}.{pf_reference.name}.final.sorted.bam"
        hs_bam_path = f"{bam_dir}/{barcode}.{hs_reference.name}.final.sorted.bam"

        # Load p.f. alignments
        print("Loading data...")
        pf_alignments_df = load_alignment_information(pf_bam_path)
        pf_alignments_df.insert(0, "species", "pf")
        pf_alignments_df.query(
            "flag != 4", inplace=True
        )  # these have been remapped to H.s.

        # Load h.s. alignments
        hs_alignments_df = load_alignment_information(hs_bam_path)
        hs_alignments_df.insert(0, "species", "hs")

        # Combine all alignments
        alignments_df = pd.concat([pf_alignments_df, hs_alignments_df])

        # Produce a read-level data frame
        print("Processing...")
        read_df = reduce_to_read_dataframe(alignments_df)
        convert_column_to_ordered_category(read_df, "primary_state", msc.primary_levels)
        convert_column_to_ordered_category(
            read_df, "secondary_state", msc.secondary_levels
        )

        # Prepare plotter
        plotter = ReadHistogramPlotter(read_df)

        # Iterate over statistics, states, and plot
        print("Plotting...")
        for state, colors in msc.color_sets.items():

            # Set states of interest, compute group sizes
            plotter.set_groups(state, colors)
            size_df = plotter.get_group_size_dataframe()
            size_df.to_csv(f"{output_dir}/table.size.{state}.csv", index=False)

            for histogram_stat in HISTOGRAM_STATS:

                # Set statistics and create histogram
                plotter.set_histogram_stats(**histogram_stat.__dict__)
                plotter.create_histogram_dataframe()

                # Plot
                plotter.plot_histogram(
                    title=barcode,
                    output_path=f"{output_dir}/plot.{histogram_stat.stat}.{state}.png",
                )

                # Write to data frame
                plotter.write_histogram_dataframe(
                    f"{output_dir}/table.bin_counts.{histogram_stat.stat}.{state}.csv"
                )
        print(f"Output directory: {output_dir}")
        print("Done.")
        print("")
    print_footer(t0)
