import re
import pandas as pd
from collections import Counter


def get_indels(pileup):
    """
    Return a list of all indels found
    in `pileup`

    params
        pileup: str
            The pileup is a string representing all
            of the bases covering a given position of
            the genome; generated using samtools mpileup.
            Indels are encoded as r"[+-]\d+[ATCG]+". e.g.
            +4AATT for an insertion of four bases.
    returns
        indels: list int
            A list of the sizes of all indels in `pileup`.
            Insertions are positive, deletions are negative.

    """

    inserts = [int(i) for i in re.findall(r"[+]\d+", pileup)]
    deletes = [int(d) for d in re.findall(r"[-]\d+", pileup)]

    return inserts + deletes


def process_mpileup(pileup, ref):
    """
    Process a read pileup from `samtools mpileup` into A,T,C,G,+,-

    Note, the solution is ugly, but re.sub([-+][0-9]+[ATCGatcg]+,
    ... fails when a variant follows directly after an indel.
    Also note that this information corresponds to an *upcoming*
    insertion or deletion; it is not to be counted as contributing
    to the present position.

    params
        pileup: str
            String giving all nucleotides mapped to a specific
            position, with characters defined by
            `samtools mpileup`.
        ref: str
            String giving the reference base at this position.
    returns
        pileup:
            String with length <= pileup. Characters given in
            `samtools mpileup` are converted to the appropriate
            A, T, C, and G values. Insertions are represented
            by a single '+', deletions by a single '-'.

    """

    # indels are initiated with [+-][0-9]+[ATCGatcg]
    for indel_size in set(re.findall("\d+", pileup)):
        indel_size = int(indel_size)
        pileup = re.sub(
            r"[+]%d[ATCGatcg]{%d}" % (indel_size, indel_size), "+", pileup
        )  # insertions
        pileup = re.sub(
            r"[-]%d[ATCGatcg]{%d}" % (indel_size, indel_size), "", pileup
        )  # deletions

    pileup = re.sub("\*", "-", pileup)  # deletions cover multiple positions
    pileup = re.sub("\$|\^.", "", pileup)  # start, or end with mapping quality
    pileup = re.sub("\.|,", ref, pileup)  # matches, . for forward strand , for reverse
    pileup = pileup.upper()

    return pileup


def create_basecall_dfs(input_pileup):
    """
    Create data frames of all the basecalls by position


    """

    nt_symbols = ["A", "T", "C", "G", "-", "+"]

    # Prepare storage
    mutation_dt = {
        "position": [],
        "ref": [],
        "coverage": [],  # Insertions not counted
        "calls": [],  # Insertions are counted
        "SNV": [],
        "error": [],
    }
    mutation_dt.update({k: [] for k in nt_symbols})

    indel_dt = {"position": [], "ref": [], "length": []}

    with open(input_pileup, "r") as fn:
        for i, line in enumerate(fn):

            # Prepare pileup
            chrom, pos, ref, coverage, pileup, _ = line.split("\t")
            processed_pileup = process_mpileup(pileup, ref)

            # Get nucleotide frequencies in pileup
            nt_frequencies = Counter(processed_pileup)

            # Compute summary statistics
            calls = sum([count for nt, count in nt_frequencies.items()])
            error = sum([count for nt, count in nt_frequencies.items() if nt != ref])
            snv = sum(
                [
                    count
                    for nt, count in nt_frequencies.items()
                    if nt != ref and not nt in ["+", "-"]
                ]
            )

            # Store
            mutation_dt["position"].append(int(pos))
            mutation_dt["ref"].append(ref)
            mutation_dt["coverage"].append(int(coverage))
            mutation_dt["calls"].append(calls)
            mutation_dt["SNV"].append(snv)
            mutation_dt["error"].append(error)
            for k in nt_symbols:
                mutation_dt[k].append(nt_frequencies[k])

            # Get indel lengths
            indels = get_indels(pileup)
            n_indels = len(indels)

            # Store
            indel_dt["position"].extend([int(pos)] * n_indels)
            indel_dt["ref"].extend([ref] * n_indels)
            indel_dt["length"].extend(indels)

    return pd.DataFrame(mutation_dt), pd.DataFrame(indel_dt)
