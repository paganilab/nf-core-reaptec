// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

process GET_WHITELIST {
    label "process_low"
    tag "$barcodes"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'genome', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    path barcodes

    output:
    path "WhiteList.txt", emit: whitelist
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    zcat $barcodes | sed -e 's/-1//g' > WhiteList.txt
    echo '0.5.0' > ${software}.version.txt
    """
}
