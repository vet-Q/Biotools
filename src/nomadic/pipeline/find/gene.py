import pysam
import pandas as pd


class Gene:
    """
    Explore codon landscape of a gene; interrogate the effect
    of SNPs on amino acid sequence

    """

    genetic_code = {
        "ATA": "I",
        "ATC": "I",
        "ATT": "I",
        "ATG": "M",
        "ACA": "T",
        "ACC": "T",
        "ACG": "T",
        "ACT": "T",
        "AAC": "N",
        "AAT": "N",
        "AAA": "K",
        "AAG": "K",
        "AGC": "S",
        "AGT": "S",
        "AGA": "R",
        "AGG": "R",
        "CTA": "L",
        "CTC": "L",
        "CTG": "L",
        "CTT": "L",
        "CCA": "P",
        "CCC": "P",
        "CCG": "P",
        "CCT": "P",
        "CAC": "H",
        "CAT": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGA": "R",
        "CGC": "R",
        "CGG": "R",
        "CGT": "R",
        "GTA": "V",
        "GTC": "V",
        "GTG": "V",
        "GTT": "V",
        "GCA": "A",
        "GCC": "A",
        "GCG": "A",
        "GCT": "A",
        "GAC": "D",
        "GAT": "D",
        "GAA": "E",
        "GAG": "E",
        "GGA": "G",
        "GGC": "G",
        "GGG": "G",
        "GGT": "G",
        "TCA": "S",
        "TCC": "S",
        "TCG": "S",
        "TCT": "S",
        "TTC": "F",
        "TTT": "F",
        "TTA": "L",
        "TTG": "L",
        "TAC": "Y",
        "TAT": "Y",
        "TAA": "_",
        "TAG": "_",
        "TGC": "C",
        "TGT": "C",
        "TGA": "_",
        "TGG": "W",
    }

    def __init__(self, gff_path, fasta_path):

        self.gff = self._load_gff(gff_path)
        self.gff = self._add_gff_fields(self.gff)
        self.fasta_path = fasta_path
        self.gene_id = None

        self.cds_table = None

        self.coding_sequence = None
        self.positions = None
        self.chrom = None
        self.strand = None
        self.L = None

        self.codon_nts = None
        self.codon_aas = None

    @staticmethod
    def _load_gff(gff_path):
        """Load .gff file into a dataframe"""

        # Count number of header lines in .gff file
        with open(gff_path, "r") as gff:
            cnt = 0
            for line in gff:
                if line.startswith("##"):
                    cnt += 1
                else:
                    break

        # Load .gff file
        gff = pd.read_csv(
            gff_path,
            sep="\t",
            skiprows=cnt,
            header=None,
            names=[
                "seqname",
                "source",
                "feature",
                "start",
                "end",
                "score",
                "strand",
                "frame",
                "attribute",
            ],
        )

        return gff

    @staticmethod
    def _add_gff_fields(gff):
        """Add name, ID, and description fields to gff data frame"""

        # Fields
        name = []
        parent = []
        ID = []
        descrip = []
        fields = [name, parent, ID, descrip]

        # Tags
        tags = ["Name=", "Parent=", "ID=", "description="]

        # Iterate over rows, extract new fields
        for ix, row in gff.iterrows():
            attributes = row["attribute"].split(";")
            for field, tag in zip(fields, tags):
                r = [a.replace(tag, "") for a in attributes if a.startswith(tag)]
                field.extend(r if r else [None])

        # Add new fields to data frame
        gff["name"] = name
        gff["parent"] = parent
        gff["ID"] = ID
        gff["descrip"] = descrip

        return gff

    @staticmethod
    def _load_haplotype_from_fasta(fasta_path, **kwargs):
        """
        Load a particular region of a fasta file

        """

        if kwargs is None:
            fasta_seq = next(pysam.FastxFile(fasta_path)).sequence
        else:
            fasta_seq = pysam.FastaFile(fasta_path).fetch(**kwargs)

        return fasta_seq

    @staticmethod
    def _reverse_compliment(seq):
        """Reverse compliment a sequence"""

        return seq.translate(seq.maketrans("ACTG-", "TGAC-"))[::-1]

    def set_gene_id(self, gene_id):
        """Sene the gene of interest by its ID"""

        self.gene_id = gene_id

        self._create_gene_coding_table()
        self._get_coding_sequence()

        return None

    def _create_gene_coding_table(self):
        """Create a table of coding regions from the .gff"""

        ixs = [
            ix
            for ix, p in self.gff["parent"].iteritems()
            if p and p.startswith(self.gene_id)
        ]

        self.cds_table = self.gff.loc[ixs].query("feature == 'CDS'")
        self.strand = str(self.cds_table.iloc[0]["strand"])
        self.chrom = str(self.cds_table.iloc[0]["seqname"])

        return None

    def _get_coding_sequence(self):
        """
        Get the coding sequence in nucleotides and amino acids
        for the target `gene_id`

        """

        # Get nucleotide sequence
        self.coding_sequence = ""
        self.positions = []
        for _, row in self.cds_table.iterrows():

            self.coding_sequence += self._load_haplotype_from_fasta(
                self.fasta_path,
                reference=row["seqname"],
                start=row["start"] - 1,
                end=row["end"],
            )

            self.positions.extend(list(range(row["start"], row["end"] + 1)))
        self.L = len(self.coding_sequence)

        # Reverse compliment, if necessaray
        if self.strand == "-":
            self.coding_sequence = self._reverse_compliment(self.coding_sequence)
            self.positions = self.positions[::-1]

        # Group into codons
        self.codon_positions = [
            self.positions[i : i + 3] for i in range(0, self.L - 2, 3)
        ]
        self.codon_nts = [
            self.coding_sequence[i : i + 3] for i in range(0, self.L - 2, 3)
        ]
        self.codon_aas = [
            self.genetic_code[nts] if nts in self.genetic_code else None
            for nts in self.codon_nts
        ]

        assert self.codon_aas[0] == "M"
        assert self.codon_aas[-1] == "_"

        return None

    def get_codon(self, number):
        """Get a codon by its number"""

        return self.codon_nts[number - 1]

    def get_aa(self, number):
        """Get an amino acid by its number"""

        return self.codon_aas[number - 1]

    def query_for_mutation(self, chrom, pos, ref, alt):
        """
        Check if a given mutation, specified by a chromosome `chrom`,
        position `pos`, reference base `ref` and alternative base `alt`
        exists within the gene.
        
        If it does, return the amino-acid level consequence of the mutation,
        e.g. K76T.
        
        If it does not, return None.
        
        Note that all variants are represented with respect to the forward strand,
        though some genes may be on the reverse strand.
        
        """

        if not self.chrom == chrom:
            return None

        if pos not in self.positions:
            return None
        
        if self.strand == "-":
            ref = self._reverse_compliment(ref)
            alt = self._reverse_compliment(alt)

        ref_found = self.coding_sequence[self.positions.index(pos)]
        assert (
            ref_found == ref
        ), f"Found mutation, but expected {ref} and found {ref_found}."

        ix = [j for j, ps in enumerate(self.codon_positions) if pos in ps][0]
        cix = self.codon_positions[ix].index(pos)

        ref_aa = self.codon_aas[ix]
        ref_nts = self.codon_nts[ix]
        alt_nts = list(ref_nts)
        alt_nts[cix] = alt
        alt_nts = "".join(alt_nts)
        alt_aa = self.genetic_code[alt_nts] if alt_nts in self.genetic_code else None

        return f"{ref_aa}{ix+1}{alt_aa}"