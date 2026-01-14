process MAP_TO_PF3D7 {
    tag "map_pf3d7"
    cpus params.threads
    publishDir "${params.outdir}/03_pf_mapping", mode: 'copy'

    // No channel input: process will read files directly from a configured demux directory
    output:
        path "mapped_pf/*.bam"
        path "mapped_pf/*.bam.bai"

    script:
        """
        bash ${projectDir}/bin/map_to_pf3d7.sh ${projectDir}/${params.demuxdir} ${projectDir}/${params.pf_fa} ${params.threads}
        """
}

workflow {
    MAP_TO_PF3D7()
}