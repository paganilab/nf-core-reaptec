// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process UNENCODED_G {
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
    tuple val(meta), path(bam)
    val(base)

    output:
    tuple val(meta), path("Softclip${base}_*.bam"), emit: clipped
    tuple val(meta), path("*.log"), emit: log
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    UnencodedG.sh ${bam} ${base} ${task.cpus} > ${software}.log 2>&1
    echo '0.0.1' > ${software}.version.txt
    """
}
