import numpy as np
import pandas as pd

from dataclasses import dataclass
from typing import List


@dataclass(frozen=False)
class Read:
    read_id: str
    seq: str
    quals: str
    length: int = None

    def __post_init__(self):
        self.length = len(self.seq)


@dataclass
class ReadInfo:
    read_id: str
    length: int
    per_gc: float
    mean_qual: float
    med_qual: float


def load_fastq_reads(fastq_path: str) -> List[Read]:
    """
    Load reads from a fastq file at `fastq_path`
    
    """
    reads = []
    with open(fastq_path, "r") as fastq:
        for line in fastq:
            if line.startswith("@"):
                name = line.strip()[1:]
                seq = fastq.readline().strip()
                assert fastq.readline().strip() == "+"
                quals = fastq.readline().strip()
                reads.append(Read(read_id=name, seq=seq, quals=quals))

    return reads


def calc_percent_gc(seq: str) -> float:
    """ Calculate GC percentage of a sequece """
    N = len(seq)
    n_gc = len([nt for nt in seq if nt in ["G", "C"]])
    return 100*n_gc/N


def convert_ascii_to_quals(ascii_quals: str) -> np.ndarray:
    """
    Convert ASCII representation of quality scores
    to numeric

    """
    return np.array([ord(q) - 33 for q in ascii_quals])


def load_fastq_read_info(fastq_path: str) -> pd.DataFrame:
    """
    Load reads from `fastq_path` and return information
    about them
    
    """
    
    # Load reads
    reads = load_fastq_reads(fastq_path=fastq_path)

    # Compute summary statistics
    results = []
    for r in reads:
        quals = convert_ascii_to_quals(r.quals)
        results.append(
            ReadInfo(
                read_id=r.read_id,
                length=r.length,
                per_gc=calc_percent_gc(r.seq),
                mean_qual=quals.mean(),
                med_qual=np.median(quals)
            )
        )
    return pd.DataFrame(results)

