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
    ),
    Target(
        name="CSP",
        chrom="Pf3D7_03_v3",
        start=221323,
        end=222516
    ),
    Target(
        name="AMA1",
        chrom="Pf3D7_11_v3",
        start=1293856,
        end=1295724
    )
]


TARGET_COLLECTION = {
    t.name: t for t in TARGETS
}