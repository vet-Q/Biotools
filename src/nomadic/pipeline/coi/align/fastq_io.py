from dataclasses import dataclass
from typing import List

@dataclass(frozen=False)
class Read:
    name: str
    seq: str
    quals: str
    seq_len: int = None

    def __post_init__(self):
        self.seq_len = len(self.seq)


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
                reads.append(Read(name=name, seq=seq, quals=quals))

    return reads