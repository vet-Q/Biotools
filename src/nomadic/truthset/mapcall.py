# Steps:
# - Load FASTAs
# - Convert to .fastq (is this necessary? does minimap2 handle .fasta)
#  - This would be done for each of the strains
# - Map with given alignment method (producing .bam)
# - Call variants with a given (or multiple?) aligment
#  - Probably just one method, and have output as sub-directory so can subsequently compare

import os
from re import S
import click
from nomadic.lib.generic import print_header, print_footer, produce_dir
from nomadic.lib.references import PlasmodiumFalciparum3D7
from nomadic.pipeline.map.mappers import MAPPER_COLLECTION
from nomadic.pipeline.calling.callers import caller_collection
from nomadic.lib.process_bams import run_samtools


# ================================================================
# Settings
#
# ================================================================


N_READS = 30
QUALITY_SCORE_CONST = 30


# ================================================================
# File conversion
#
# ================================================================


def convert_fasta_to_fastq(fasta_path, n_reads, quality_score_const, fastq_path):
    """
    Convert a FASTA file to a FASTQ file, with `n_reads`
    per sequence and a constant quality score `quality_score_const`

    """
    # Phred score as ASCII
    quality_ascii = chr(quality_score_const + 33)

    # Open FASTQ
    with open(fastq_path, "w") as fastq:

        # Open FASTA
        with open(fasta_path, "r") as fasta:

            # Iterate over FASTQ, producing reads in FASTQ
            for header in fasta:

                # Extract header and sequence
                assert header.startswith(">")
                strain = header.strip()[1:]
                sequence = fasta.readline().strip()

                # Reverse complement sequence
                tt = sequence.maketrans({"A": "T", "T": "A", "C": "G", "G": "C"})
                rc_sequence = sequence.translate(tt)[::-1]

                # Template for a single FASTQ read
                fastq_read = "@{strain} strand={strand} read={n_rep:05d}\n"
                fastq_read += "{sequence}\n"
                fastq_read += "+\n"
                fastq_read += "{quals}\n"

                # Write for forward strand
                fastq.write(
                    "".join(
                        [
                            fastq_read.format(
                                strain=strain,
                                strand=strand,
                                n_rep=ix,
                                sequence=strand_sequence,
                                quals=quality_ascii * len(sequence),
                            )
                            for strand, strand_sequence in {
                                "F": sequence,
                                "R": rc_sequence,
                            }.items()
                            for ix in range(n_reads)
                        ]
                    )
                )


# ================================================================
# Main script
#
# ================================================================


@click.command(short_help="Create a truthset via mapping and calling.")
@click.option(
    "-f",
    "--fasta_dir",
    type=click.Path(exists=True),
    default="resources/truthsets/fastas",
    help="Path to directory containing .fasta files.",
)
@click.option(
    "-m",
    "--map_method",
    type=click.Choice(MAPPER_COLLECTION),
    default="minimap2",
    help="Algorithm used to map reads.",
)
@click.option(
    "-c",
    "--call_method",
    type=click.Choice(caller_collection),
    default="bcftools",
    help="Algorithm used to call variants.",
)
def mapcall(fasta_dir, map_method, call_method):
    """
    Create a truth set of variants by mapping FASTA files from PacBio assemblies
    to the 3D7 reference, and calling variants

    """
    # Parse CLI
    script_descrip = (
        "TRUTHSET: Create a truth variant set using read mapping and variant calling"
    )
    script_name = "mapcall"
    t0 = print_header(script_descrip)
    truthset_dir = "resources/truthsets"
    output_dir = produce_dir(truthset_dir, script_name)
    print("Parsing inputs...")
    print(f"  Mapping method: {map_method}")
    print(f"  Calling method: {call_method}")
    print("Done.")
    print("")

    # Collect FASTA files from `fasta_dir`
    print("Preprocessing FASTA files...")
    fastas = [
        f"{fasta_dir}/{f}"
        for f in os.listdir(fasta_dir)
        if f.endswith(".fasta") or f.endswith(".fa")
    ]
    print(f"  FASTA directory: {fasta_dir}")
    print(f"  Found {len(fastas)} FASTA files.")

    # Populate a dictionary of per-strain sequences
    print("Creating per-strain fasta files...")
    strain_fasta = {}
    for fasta_path in fastas:
        # Store sequences in `strain_fasta`
        with open(fasta_path, "r") as f:
            for header in f:
                assert header.startswith(">")
                strain = header.split("|")[0].strip()[1:]
                gene_id = header.split("|")[1].strip()
                sequence = f.readline().strip()
                if strain in strain_fasta:
                    strain_fasta[strain][gene_id] = sequence
                else:
                    strain_fasta[strain] = {gene_id: sequence}

    # Sanity checks...
    print(f"  Found{len(strain_fasta)} strains: {', '.join(strain_fasta)}")
    n_seqs = [len(v) for k, v in strain_fasta.items()]
    assert min(n_seqs) == max(
        n_seqs
    ), "Not all strains have the same number of associated sequences."
    print(f"  Each strain has {max(n_seqs)} sequences.")

    # Write fasta file for each strain
    output_fasta_dir = produce_dir(output_dir, "fastas")
    print(f"  Writing output fastas to: {output_fasta_dir}")
    strain_fasta_paths = []
    for strain, gene_sequences in strain_fasta.items():

        # Define an output fasta for the strain
        strain_fasta_path = f"{output_fasta_dir}/{strain}.fasta"
        strain_fasta_paths.append(strain_fasta_path)

        # Write file
        with open(strain_fasta_path, "w") as f:
            for gene_id, sequence in gene_sequences.items():
                f.write(f">{gene_id}\n")
                f.write(f"{sequence}\n")
    print("Done.")
    print("")

    #  Create FASTQs
    fastq_output_dir = produce_dir(output_dir, "fastqs")
    fastq_output_paths = []
    for strain_fasta_path in strain_fasta_paths:
        fastq_path = f"{fastq_output_dir}/{os.path.basename(strain_fasta_path).replace('.fasta', '.fastq')}"
        fastq_output_paths.append(fastq_path)
        convert_fasta_to_fastq(
            fasta_path=strain_fasta_path,
            n_reads=N_READS,
            quality_score_const=QUALITY_SCORE_CONST,
            fastq_path=fastq_path,
        )

    # Mapping
    reference = PlasmodiumFalciparum3D7()
    print(f"Mapping FASTQ files to {reference.name}...")
    mapper = MAPPER_COLLECTION[map_method](reference)
    bam_output_dir = produce_dir(output_dir, "bams", map_method)
    output_bams = []
    for fastq_path in fastq_output_paths:

        # FASTQ
        print(f"Input fastq: {fastq_path}")

        # Map
        output_bam = (
            f"{bam_output_dir}/{os.path.basename(fastq_path).replace('.fastq', '.sam')}"
        )
        print("Mapping...")
        mapper.map_from_fastqs(fastq_path=fastq_path)
        mapper.run(output_bam)

        # Sort
        print("  Sorting...")
        output_sorted_bam = output_bam.replace(".sam", ".sorted.sam")
        run_samtools(utility="sort", args=f"{output_bam} -o {output_sorted_bam}")
        print("Done.")
        print("")

        # Store
        output_bams.append(output_sorted_bam)
    print("Done.")
    print("")

    # Calling
    print("Calling variants...")
    caller = caller_collection[call_method]
    output_vcf_dir = produce_dir(output_dir, "vcfs", call_method)
    for bam_path in output_bams:

        # Define output VCF file
        strain, _, _ = os.path.basename(bam_path).split(".")
        vcf_path = f"{output_vcf_dir}/{strain}.{map_method}.vcf"
        print(f"Input BAM: {bam_path}")
        print(f"Output VCF: {vcf_path}")

        # Call variants
        caller.set_files(bam_path=bam_path, vcf_path=vcf_path)
        caller.set_arguments(fasta_path=reference.fasta_path)
        caller.call_variants(sample_name=strain)
        print("Done.")
    print("Done.")
    print("")

    print_footer(t0)
