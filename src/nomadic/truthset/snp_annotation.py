import pandas as pd
from collections import namedtuple
from src.nomadic.lib.process_gffs import load_gff, add_gff_fields


class SNPAnnotator:
    def __init__(self, gff):
        """ 
        Create Annotation Data Frames based on Gene IDs;
        Use them to annotate a set of SNPs based on their
        position
        
        """
        
        for c in ["name", "ID"]:
            assert c in gff.columns
            
        self.gff = gff
        self.gene_ids = None
        
        self.genome_df = None    
        self.amplicon_df = None
        
        self.cds_sizes = None
        self.amplicon_sizes = None
    
    
    @classmethod
    def from_gff(cls, gff_path):
        """ Instantiate the class from a gff """
        
        gff = load_gff(gff_path)
        gff = add_gff_fields(gff)
        gff.rename({"Name": "name", "Parent": "parent"}, axis=1, inplace=True)
        gff["start"] = gff["start"].astype("int")
        gff["end"] = gff["end"].astype("int")
        
        return cls(gff)
        
    
    def provide_gene_ids(self, gene_ids):
        """ Provide a list of gene IDs of interest """
        
        self.gene_ids = gene_ids
            
        return None
        
        
    def create_genome_table(self):
        """ Create a table of possible annotations """

        dfs = []
        for gene_id in self.gene_ids:
            idxs = [j for j, r in enumerate(self.gff["ID"]) if gene_id in r]
            if idxs:
                dfs.append(self.gff.iloc[idxs])
            else:
                raise Warning(f"No .gff entries found for {gene_id}.")
        df = pd.concat(dfs)

        # Filter to key features
        keep_features = ["CDS", "five_prime_UTR", "three_prime_UTR"]
        self.genome_df = df.query("feature in @keep_features")
        self._get_annotation_sizes()

        return None
    
    
    def create_amplicon_table(self, csv_path):
        """ 
        Load data about multiplex being used, 
        convert to annotation table

        """

        # Load
        df = pd.read_csv(csv_path)

        # Build amplicons
        Amplicon = namedtuple("Amplicon", 
                              ["ID", "name", "seqname", "start", "end"])
        amplicons = []
        for gene_id in self.gene_ids:
            adf = df.query("target == @gene_id")
            if not adf.shape[0] == 2:
                print(f"Failed to find amplicon for {gene_id}. Skipping.")
                continue
            chrom_num = int(gene_id.split("_")[1][:2])
            amplicon = Amplicon(
                gene_id, 
                adf.iloc[0]["gene_name"],
                f"Pf3D7_{chrom_num:02d}_v3",
                adf.iloc[0]["position"] - 1,
                adf.iloc[1]["position"]
            )
            assert amplicon.start < amplicon.end
            amplicons.append(amplicon)

        # Save as dataframe
        self.amplicon_df = pd.DataFrame(amplicons)
        #self._get_annotation_sizes()

        return None
    

    def _get_annotation_sizes(self):
        """ Get sizes of annotations currently loaded """ 

        if self.genome_df is not None:

            cds_df = self.genome_df.query("feature == 'CDS'")
            self.cds_sizes = {}

            for gene_id, gdf in cds_df.groupby("parent"):

                if "." in gene_id:
                    gene_id = gene_id.split(".")[0]

                size = (gdf["end"] - gdf["start"]).sum()
                self.cds_sizes[gene_id] = size

        if self.amplicon_df is not None:
            self.amplicon_sizes = {
                r["ID"]: int(r["end"] - r["start"])
                for _, r in self.amplicon_df.iterrows()
            }

        return None

    
    def _find_annotation(self, df, chrom, pos, annot_field="feature"):
        """ Find annotations in a dataframe `df` """
        
        search = df.query(
            "seqid == @chrom and start <= @pos <= end"
        )
        
        if search.shape[0] == 1:
            annot = search.iloc[0][annot_field]
        elif search.shape[0] == 0:
            annot = "None"
        else:
            raise Warning("More than one annotation found for this SNP!")
            
        return annot
        
    
    def annotate_snp_by_genome(self, chrom, pos, annot_field="feature"):
        """ 
        Annotate where a SNP is within the genome
        based on its chromosome and position
        
        """
        
        return self._find_annotation(self.genome_df, 
                                     chrom, pos, annot_field)
    
    
    def annotate_snp_by_amplicon(self, chrom, pos):
        """ 
        Annotate whether a SNP is inside an amplicion 
        based on its chromosome and position
        
        """
        
        return self._find_annotation(self.amplicon_df, 
                                     chrom, pos, 
                                     annot_field="name")