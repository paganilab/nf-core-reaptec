// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BAM_TO_CTSS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    // TODO: verify the right container and the requirements
    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    path(index) // STAR index path
    path(promoters_bed)
    path(enhancers_bed)
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*_mq20.ctss.bed.gz")                          , emit: ctss
    tuple val(meta), path("*_promoter.count.txt")                        , emit: promoters
    tuple val(meta), path("*_enhancer.count.txt")                        , emit: enhancers
    tuple val(meta), path("*.fwd.bw"), path("*.rev.bw"), path("*.all.bw"), emit: bigwigs
    tuple val(meta), path("*.log"), emit: log
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    BAMtoCTSS.sh ${index}/chrNameLength.txt $promoters_bed $enhancers_bed $bam ${task.cpus} > ${software}.log 2>&1
    echo '0.0.1' > ${software}.version.txt
    """
}
