// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CALL_ENHANCERS {
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
    tuple val(meta), path(ctss)
    path(mask)
    val(prefix)

    output:
    tuple val(meta), path("output"), emit: enhancers
    tuple val(meta), path("*.log"), emit: log
    path "*.version.txt", emit: version

    script:
    def software = 'call_enhancers' //getSoftwareName(task.process)
    """
    gunzip -c ${ctss} > ctss.bed
    echo "ctss.bed" > ctss.txt
    fixed_bidir_enhancers_10bp -s ${prefix}_ -m $mask -t tmp -f ctss.txt -o output > ${software}.log 2>&1
    echo '0.0.1' > ${software}.version.txt
    """
}
