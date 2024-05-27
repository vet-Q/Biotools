import click
from nomadic.pipeline.calling.downsample import BamDownSampler
from nomadic.lib.parsing import build_parameter_dict
from nomadic.pipeline.cli import experiment_options, barcode_option


N_READS = 20_000


@click.command(short_help="Downsample BAM files.")
@experiment_options
@barcode_option
@click.option(
    "-n",
    "--n_reads",
    type=int,
    default=N_READS,
    show_default=True,
    help="Downsample to this many reads."
)
def main(expt_dir: str, config: str, barcode: str = None, n_reads: int = N_READS):
    """
    Downsample all of the BAM files in an experiment to a prespecified depth
    
    """

    # PARSE INPUTS
    params = build_parameter_dict(expt_dir, config, barcode)

    # Focus on a single barcode, if specified
    if "focus_barcode" in params:
        params["barcodes"] = [params["focus_barcode"]]

    # ITERATE OVER BAMS
    for barcode in params["barcodes"]:
        bam_dir = f"{params['barcodes_dir']}/{barcode}/bams"
        input_bam = f"{bam_dir}/{barcode}.Pf3D7.final.sorted.bam"
        output_bam = input_bam.replace(".bam", ".downsampled.bam")

        # Initialise downsampler
        downsampler = BamDownSampler(input_bam)
        if downsampler.n_reads_total < n_reads:
            print(f"Not enough reads for {barcode} to downsample.")
            print(f"  No. reads total: {downsampler.n_reads_total}")
            print(f"  No. downsampling: {n_reads}")
            continue

        print("Downsampling...")
        downsampler.create_downsample(n_reads=n_reads,
                                      downsampled_bam_path=output_bam
                                      )
        print("Done.")
        print("")


if __name__ == "__main__":
    main()