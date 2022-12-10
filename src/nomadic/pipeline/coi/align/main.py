import datetime
import numpy as np
import pandas as pd

from nomadic.lib.generic import produce_dir, print_header, print_footer
from nomadic.lib.parsing import build_parameter_dict
from nomadic.lib.process_fastqs import load_fastq_reads
from ..trim.targets import TARGET_COLLECTION
from .aligners import ALIGNER_COLLECTION


def main(expt_dir, config, barcode, target_gene, max_reads, algorithm):

    # PARSE INPUTS
    script_descrip = "NOMADIC: Map reads from a target gene to a panel of P.f. strains."
    t0 = print_header(script_descrip)
    script_dir = "coi"
    params = build_parameter_dict(expt_dir, config, barcode)

    target = TARGET_COLLECTION[target_gene]
    aligner = ALIGNER_COLLECTION[algorithm]()
    print("User inputs:")
    print(f"  Target: {target.name}")
    print(f"  Chrom: {target.chrom}")
    print(f"  Start: {target.start}")
    print(f"  End: {target.end}")
    print(f"  Alignment algorithm: {algorithm}")
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

        bt0 = datetime.datetime.now().replace(microsecond=0)

        # DIRECTORIES
        coi_dir = produce_dir(params["barcodes_dir"], barcode, script_dir)
        fastq_path = f"{coi_dir}/fastq_clipped/reads.target.{target_gene}.clipped.fastq"
        output_dir = produce_dir(coi_dir, "align")

        # Load FASTQ
        print("Loading reads from FASTQ...")
        reads = load_fastq_reads(fastq_path)
        n_reads = len(reads)
        print(f" Found {n_reads} reads...")


        if n_reads > max_reads:
            print(f"  Exceeds maximum of {max_reads}!")
            print(f"  Reducing to first {max_reads}.")
            reads = reads[:max_reads]
            n_reads = max_reads

        print("Performing pairwise alignments...")
        scores = np.zeros((n_reads, n_reads))
        for i in range(n_reads):
            for j in range(i, n_reads):
                aligner.set_sequences(
                    x=reads[i].seq, 
                    y=reads[j].seq
                )
                aligner.set_scoring_model()
                aligner.align()
                scores[i, j] = aligner.score
                scores[j, i] = aligner.score
                
            if i % 100 == 0:
                print(f"Completed first {i} reads...")

        read_names = [r.read_id for r in reads]
        score_df = pd.DataFrame(
            scores,
            index=read_names,
            columns=read_names
        )
        score_df.to_csv(f"{output_dir}/pairwise_scores.{target_gene}.csv")

        bt1 = datetime.datetime.now().replace(microsecond=0)
        print("Time Elapsed: %s" % (bt1 - bt0))
    
    print_footer(t0)