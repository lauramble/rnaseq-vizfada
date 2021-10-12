// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CAT_FASTQ {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'merged_fastq', meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }
    
    afterScript "if [ ${params.keep_fastq} == 'false' ];\
                 then \
                    ls -l *.fastq.gz & \
                    ls -l *.fastq.gz | \
                    grep -- '->' | \
                    sed -e's/.*-> //' | \
                    xargs rm;\
                 fi;"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.merged.fastq.gz"), emit: reads

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def readList = reads.collect{ it.toString() }
    if (meta.single_end) {
        if (readList.size > 1) {
            """
            cat ${readList.sort().join(' ')} > ${prefix}.merged.fastq.gz
            """
        }
    } else {
        if (readList.size > 2) {
            def read1 = []
            def read2 = []
            readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }
            """
            cat ${read1.sort().join(' ')} > ${prefix}_1.merged.fastq.gz
            cat ${read2.sort().join(' ')} > ${prefix}_2.merged.fastq.gz
            """
        }
    }
}
