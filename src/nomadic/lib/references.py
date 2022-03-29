import os
import urllib.request
from abc import ABC, abstractmethod


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
    release = 52

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
    source_url += "000/001/405/GCA_000001405.15_GRCh38/GO_TO_CURRENT_VERSION/"

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
            print("Downloading...")
            self.produce_dir(self.ref.fasta_path)
            urllib.request.urlretrieve(
                url=self.ref.fasta_url, filename=self.ref.fasta_path
            )
            print("Done.")
            print("")
        else:
            print("Already downloaded.")

    def download_gff(self):
        if self.ref.gff_path and not self.exists_locally(self.ref.gff_path):
            print("Downloading...")
            self.produce_dir(self.ref.gff_path)
            urllib.request.urlretrieve(url=self.ref.gff_url, filename=self.ref.gff_path)
            print("Done.")
        else:
            print("Already downloaded.")

            
if __name__ == "__main__":
    references = [
        PlasmodiumFalciparum3D7(),
        PlasmodiumFalciparumDd2(),
        PlasmodiumOvale(),
        PlasmodiumMalariae(),
        HomoSapiens()
    ]
    downloader = ReferenceDownloader()
    for r in references:
        print(f"Downloading: {r.name}")
        downloader.set_reference(r)
        downloader.download_fasta()
        downloader.download_gff()
        print("Done.")
        print("")
            