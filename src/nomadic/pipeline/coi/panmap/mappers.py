from nomadic.pipeline.map.mappers import MappingAlgorithm


class Minimap2PAF(MappingAlgorithm):
    """
    Map long reads with `minimap2`

    """

    def _define_mapping_command(self, output_bam):
        """
        Run minimap2, compress result to .bam file, and sort

        """

        # TO OUTPUT AS PAF
        self.map_cmd = "minimap2"
        self.map_cmd += f" -x map-ont {self.reference.fasta_path} {self.input_fastqs}"
        self.map_cmd += f" > {output_bam}"

