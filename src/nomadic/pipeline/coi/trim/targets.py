from dataclasses import dataclass


@dataclass
class Target:
    name: str
    chrom: str
    start: int
    end: int


TARGETS = [
    Target(
        name="MSP2",
        chrom="Pf3D7_02_v3",
        start=273689,
        end=274507
    )
]


TARGET_COLLECTION = {
    t.name: t for t in TARGETS
}