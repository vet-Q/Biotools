import os
import urllib.request
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List
from nomadic.lib.process_gffs import load_gff
from nomadic.lib.process_fastas import load_fasta_as_dict, write_fasta_from_dict


# ================================================================
# Parent classes for reference sequences
#
# ================================================================


class Reference(ABC):
    """
    Basic interface for reference sequences used for
    mapping, variant calling, &c.

    """

    def __init__(self):
        self.name = None
        self.fasta_url = None
        self.gff_url = None
        self.fasta_path = None
        self.gff_path = None

    @property
    def gff_standard_path(self):
        if self.gff_path is None:
            return None
        return self.gff_path.replace(".gff", ".standardised.gff")

    @abstractmethod
    def set_fasta(self):
        pass

    @abstractmethod
    def set_gff(self):
        pass


class PlasmoDB(Reference):
    """
    Encapsulate reference sequence downloads from PlasmoDB

    """

    source = "plasmodb"
    source_url = "https://plasmodb.org/common/downloads"
    release = 67

    def __init__(self, species, strain):
        self.species = species
        self.strain = strain
        self.data_url = f"{self.source_url}/release-{self.release}/{species}{strain}"
        self.set_fasta()
        self.set_gff()

    def set_fasta(self):
        """Set .fasta file download URL and local path"""
        fasta_fn = f"PlasmoDB-{self.release}_{self.species}{self.strain}_Genome.fasta"
        self.fasta_url = f"{self.data_url}/fasta/data/{fasta_fn}"
        self.fasta_path = f"resources/{self.source}/{self.release}/{fasta_fn}"

    def set_gff(self):
        """Set .gff file download URL and local path"""

        gff_fn = f"PlasmoDB-{self.release}_{self.species}{self.strain}.gff"
        self.gff_url = f"{self.data_url}/gff/data/{gff_fn}"
        self.gff_path = f"resources/{self.source}/{self.release}/{gff_fn}"


class ENA(Reference):
    """
    Encapsulate reference sequence downloads from the European Nucleotide Archive

    """

    source = "ena"
    source_url = "ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/flr"

    def __init__(self, wgs_id):
        self.wgs_id = wgs_id
        self.set_fasta()
        self.set_gff()

    def set_fasta(self):
        fasta_fn = f"{self.wgs_id}.fasta.gz"
        self.fasta_url = f"{self.source_url}/{fasta_fn}"
        self.fasta_path = f"resources/{self.source}/{fasta_fn}"

    def set_gff(self):
        self.gff_url = None
        self.gff_path = None


# ===============================================================
# Classes for specific reference sequences
#
# ================================================================


class PlasmodiumFalciparum3D7(PlasmoDB):
    def __init__(self):
        self.name = "Pf3D7"
        super().__init__(species="Pfalciparum", strain="3D7")


class PlasmodiumFalciparumDd2(PlasmoDB):
    def __init__(self):
        self.name = "PfDd2"
        super().__init__(species="Pfalciparum", strain="Dd2")


class PlasmodiumFalciparumHB3(PlasmoDB):
    def __init__(self):
        self.name = "PfHB3"
        super().__init__(species="Pfalciparum", strain="HB3")


class PlasmodiumFalciparumGB4(PlasmoDB):
    def __init__(self):
        self.name = "PfGB4"
        super().__init__(species="Pfalciparum", strain="GB4")


class PlasmodiumOvale(ENA):
    def __init__(self):
        self.name = "Po"
        super().__init__(wgs_id="FLRI01")


class PlasmodiumMalariae(ENA):
    def __init__(self):
        self.name = "Pm"
        super().__init__(wgs_id="FLRK01")


class HomoSapiens(Reference):
    """
    Download the Homo Sapiens reference genome from
    Ensembl

    """

    source = "ensembl"
    source_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/"
    source_url += "000/001/405/GCA_000001405.29_GRCh38.p14/"

    def __init__(self):
        self.name = "Hs"
        self.set_fasta()
        self.set_gff()

    def set_fasta(self):
        fasta_fn = "GCA_000001405.29_GRCh38.p14_genomic.fna.gz"
        self.fasta_url = f"{self.source_url}/{fasta_fn}"
        self.fasta_path = f"resources/{self.source}/{fasta_fn}"

    def set_gff(self):
        gff_fn = "GCA_000001405.29_GRCh38.p14_genomic.gff.gz"
        self.gff_url = f"{self.source_url}/{gff_fn}"
        self.gff_path = f"resources/{self.source}/{gff_fn}"


# ================================================================
# Downloader for specific reference sequences
#
# ================================================================


class ReferenceDownloader:
    def __init__(self):
        self.ref = None

    def set_reference(self, reference):
        self.ref = reference

    @staticmethod
    def exists_locally(file_path):
        return os.path.isfile(file_path)

    @staticmethod
    def produce_dir(file_path):
        file_dir = os.path.dirname(file_path)
        if not os.path.isdir(file_dir):
            os.makedirs(file_dir)

    def download_fasta(self):
        if self.ref.fasta_path and not self.exists_locally(self.ref.fasta_path):
            print("Downloading FASTA...")
            self.produce_dir(self.ref.fasta_path)
            urllib.request.urlretrieve(
                url=self.ref.fasta_url, filename=self.ref.fasta_path
            )
            print("Done.")
            print("")
        else:
            print("Already downloaded.")

    def download_gff(self, standardise: bool = True):
        if self.ref.gff_path and not self.exists_locally(self.ref.gff_path):
            print("Downloading GFF...")
            self.produce_dir(self.ref.gff_path)
            urllib.request.urlretrieve(url=self.ref.gff_url, filename=self.ref.gff_path)
            print("Done.")
        else:
            print("Already downloaded.")
        if standardise:
            print("Standardising GFF...")
            self._standardise_gff(self.ref.gff_path)

    def _standardise_gff(self, gff_path: str):
        """
        Try to standardise the GFF file into GFF3 format

        """
        
        # Settings
        KEEP_FIELDS = ["protein_coding_gene", "mRNA", "exon", "CDS"]
        to_gff3 = {
            "protein_coding_gene": "gene",
            "mRNA": "transcript"
        }

        # Standardise
        gff_df = load_gff(gff_path)
        gff_df.query("feature in @KEEP_FIELDS", inplace=True)
        gff_df["feature"] = [
            to_gff3[f] if f in to_gff3 else f
            for f in gff_df["feature"]
        ]

        # Write to 'standardised' path
        gff_df.to_csv(self.ref.gff_standard_path, sep="\t", index=False, header=False)


# ================================================================
# Convert DHPS to wild-type in 3D7
#
# ================================================================

@dataclass
class NucleotideChange:
    name: str
    chrom: str
    position: int
    before: str
    after: str

# I want to revert G>C at this position so that 3D7 carries
# an alanine (WT) rather than a glycine (Sulfadoxine R associated)
dhps = NucleotideChange(
    name="dhps-A437G", 
    chrom="Pf3D7_08_v3", 
    position=549685, 
    before="G", 
    after="C"
)

def update_reference_genome(fasta_path: str, mutations: List[NucleotideChange]) -> None:
    """
    Update a reference genome by reverting mutations; do this IN PLACE
    
    """
    dt = load_fasta_as_dict(fasta_path)

    # Iterate over mutations
    for mutation in mutations:
        print(f"Correcting: {mutation}...")
        chrom_str = ">" + mutation.chrom

        found_chrom = False
        for key, seq in dt.items():
            if key.startswith(chrom_str):
                found_chrom = True
                break
                

        if not found_chrom:
            sep = '\n'
            raise ValueError(f"Cannot find chromosome {mutation.chrom} in: {sep.join(dt.keys())}!")
        
        # Replace
        L = len(seq)
        dt[key] = seq[:(mutation.position - 1)] + mutation.after + seq[mutation.position:]
        assert len(dt[key]) == L, f"Something wrong with replacement, sequence length changed {L} != {len(dt[key])}."
        print("Done.")

    print("Overwriting FASTA file...")
    write_fasta_from_dict(dt, fasta_path)


# ================================================================
# Reference collection
#
# ================================================================


references = [
    PlasmodiumFalciparum3D7(),
    PlasmodiumFalciparumDd2(),
    PlasmodiumFalciparumHB3(),
    PlasmodiumFalciparumGB4(),
    #PlasmodiumMalariae(),
    #PlasmodiumOvale(),
    #HomoSapiens()
]
reference_collection = {
    ref.name: ref for ref in references
}


# ================================================================
# P.falciparum strain colors
#
# ================================================================


# PALETTE
blue = (66/255, 133/255, 244/255)
green = (52/255, 168/255, 83/255)
yellow = (251/255, 188/255, 5/255)
red = (234/255, 67/255, 53/255)
PF_REF_PALETTE = [blue, green, yellow, red]


# ================================================================
# Download references, if run as a script
#
# ================================================================


if __name__ == "__main__":
    downloader = ReferenceDownloader()
    for r in references:
        print(f"Downloading: {r.name}")
        downloader.set_reference(r)
        downloader.download_fasta()
        downloader.download_gff()
        print("Done.")
        print("")

        if r.name == 'Pf3D7':
            update_reference_genome(r.fasta_path, [dhps]) # just fix DHPS issue

            