import pysam
import pandas as pd
import numpy as np
from dataclasses import dataclass, fields


# ================================================================
# Load information about alignments from a .bam file
#
# ================================================================


def calc_percent_gc(seq):
    """Calculate GC percentage of a sequence"""

    if not seq:
        return None

    N = len(seq)
    n_gc = len([nt for nt in seq if nt in ["G", "C"]])

    return 100 * n_gc / N


def load_alignment_information(input_bam: str) -> pd.DataFrame:
    """
    Load information about every alignment from an input bam
    file `input_bam`

    """

    @dataclass
    class AlignmentSummary:
        """
        Simple summary of a single alignment
        from a .bam file

        """

        read_id: str
        mapq: int
        flag: int
        query_length: int
        query_alignment_length: int
        mean_qscore: float
        per_gc: float

        @classmethod
        def from_pysam_aligned_segment(cls, pysam_segment: pysam.AlignedSegment):
            """
            Instantiate from a pysam aligned segment

            """
            # Compute mean quality score
            qscores = np.array(pysam_segment.query_qualities)
            mean_qscore = qscores.mean() if qscores.shape else None

            return cls(
                read_id=pysam_segment.query_name,
                mapq=pysam_segment.mapq,
                flag=pysam_segment.flag,
                query_length=pysam_segment.query_length,
                query_alignment_length=pysam_segment.query_alignment_length,
                mean_qscore=mean_qscore,
                per_gc=calc_percent_gc(pysam_segment.query),
            )

    # Iterate over aligned segments in .bam, store results
    with pysam.AlignmentFile(input_bam, "r") as bam:
        results = [
            AlignmentSummary.from_pysam_aligned_segment(segment) for segment in bam
        ]

    # We need to make sure, in the case of no reads, we still produce a dataframe
    # with correct column names and types; or types get coerced later on
    if not results:
        column_types = {f.name: f.type for f in fields(AlignmentSummary)}
        empty_df = pd.DataFrame(results, columns=column_types).astype(column_types)
        return empty_df

    return pd.DataFrame(results)
