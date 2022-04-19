import pysam
import pandas as pd
from dataclasses import dataclass
from collections import namedtuple

@dataclass
class VCFRecord:
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: float
    gt: str

    @classmethod
    def from_pysam_variantfile_record(cls, record):
        pass







def load_vcf_as_df(vcf_path, only_snps=True):
    """ Load a VCF file and return as pandas DataFrame """
    
    # Define information to extract
    # - I want to add genotyping information here
    # - Probably also (maybe) want allelic depth information
    Record = namedtuple("VCFRecord", 
                        ["chrom", "pos", "ref", "alt", "qual"])
    
    # Extract
    with pysam.VariantFile(vcf_path) as vcf:
        
        records = []
        for r in vcf:
            
            if only_snps:
                
                if "INDEL" in r.info:
                    if r.info["INDEL"]:
                        continue
                        
                if len(r.ref) > 1:
                    continue

            # Store
            records.append(Record(r.chrom, int(r.pos), r.ref.upper(), r.alts[0].upper(), r.qual))
        
    return pd.DataFrame(records)




class VCFChecker:
    def __init__(self, vcf_df, mutation_df, Gene):
        """
        Check for mutations in a VCF
        
        """
        self.vcf_df = vcf_df
        self.mutation_df = mutation_df
        self.gene = Gene