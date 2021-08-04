/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowReaptec.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input,
			  params.multiqc_config,
			  params.fasta,
			  params.barcodes ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.barcodes) { ch_barcodes = file(params.barcodes) } else { exit 1, 'Input barcodes not specified!' }
/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )


/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''


//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC  } from '../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )
include { UMITOOLS_EXTRACT } from '../modules/nf-core/modules/umitools/extract/main' addParams( options: modules['umitools_extract'])
include { CUTADAPT } from '../modules/nf-core/modules/cutadapt/main' addParams( options: modules['cutadapt'] )
include { GET_WHITELIST } from '../modules/local/get_whitelist' addParams( options: [:] )
include { STAR_GENOMEGENERATE } from '../modules/nf-core/modules/star/genomegenerate/main' addParams( options: modules['star_genomegenerate'])
include { STAR_ALIGN } from '../modules/nf-core/modules/star/align/main' addParams( options: modules['star_align'])
include { SAMTOOLS_INDEX } from '../modules/nf-core/modules/samtools/index/main' addParams( options: [:])
include { UMITOOLS_DEDUP } from '../modules/nf-core/modules/umitools/dedup/main' addParams( options: modules['umitools_dedup'])
include { UNENCODED_G } from '../modules/local/unencoded_g' addParams( options: [:] )
/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow REAPTEC {

    ch_software_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))


    /* 1) Get whitelist.txt (Cell barcode list) for umi-tools
Although we can get the whitelist using umi-tools, I use the barcode list provided by CellRanger 
to proceed with the analysis using same cells that CellRanger recognizes as appropriate cells.
(This time, I used cellranger-5.0.1)
ex.
zcat ./filtered_feature_bc_matrix/barcodes.tsv.gz | sed -e 's/-1//g' > Test_Whitelist.txt
     */

    GET_WHITELIST (
	ch_barcodes
    )
    ch_software_versions = ch_software_versions.mix(GET_WHITELIST.out.version.ifEmpty(null))
	
    /* 2) Extract Cell Barcode and UMI with umi-tools from fastq files
##UMI-tools: 1.0.1 (I use this version)
     */

    UMITOOLS_EXTRACT (
	INPUT_CHECK.out.reads,
	GET_WHITELIST.out.whitelist
    )
    ch_software_versions = ch_software_versions.mix(UMITOOLS_EXTRACT.out.version.ifEmpty(null))

    // tranform the data from paired reads to single reads keeping only R1
    UMITOOLS_EXTRACT
	.out
	.reads
	.map {
	    meta, fastq ->
	    def paired_to_single_meta = [:]
	    paired_to_single_meta.id = meta.id
	    paired_to_single_meta.single_end = true
	    [paired_to_single_meta, [fastq[0]]] }
	.set { ch_just_r1 }

//3) Trim the TSO sequence (13 nt) with cutadapt from Read1 processed by umi-tools.
// #3) Trim the TSO sequence (13 nt) with cutadapt from Read1 processed by umi-tools.
// cutadapt -u 13 -o /local/home/ubuntu/DATA/ForIFOM/Process_Test_S1_L001_R1_001_trim13.fastq.gz /local/home/ubuntu/DATA/ForIFOM/Process_Test_S1_L001_R1_001.fastq.gz -j 20
    
    CUTADAPT (
	ch_just_r1
    )
    ch_software_versions = ch_software_versions.mix(CUTADAPT.out.version.ifEmpty(null))


// #4) Map read1s only using STAR
// STAR --runThreadN 12 --genomeDir /local/home/ubuntu/Ref/homSap_hg38_GENCODEv34/indexes_STAR-2.6.0c_homSap_GRCh38.p13_GENCODEv34_cmp_chr/ \
// --readFilesIn /local/home/ubuntu/DATA/ForIFOM/Process_Test_S1_L001_R1_001_trim13.fastq.gz \
// --readFilesCommand zcat \
// --outFilterMultimapNmax 1 --outTmpDir /local/home/ubuntu/DATA/ForIFOM/STAR_Tmp/ForIFOM_ \
    // --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /local/home/ubuntu/DATA/ForIFOM/STAR_results/ForIFOM_

    STAR_GENOMEGENERATE (
	params.fasta,
	params.gtf
    )
    ch_software_versions = ch_software_versions.mix(STAR_GENOMEGENERATE.out.version.ifEmpty(null))


    STAR_ALIGN (
	CUTADAPT.out.reads,
	STAR_GENOMEGENERATE.out.index,
	params.gtf
    )
    ch_software_versions = ch_software_versions.mix(STAR_ALIGN.out.version.ifEmpty(null))


    SAMTOOLS_INDEX (
	STAR_ALIGN.out.bam
    )
    ch_software_versions = ch_software_versions.mix(SAMTOOLS_INDEX.out.version.ifEmpty(null))

    UMITOOLS_DEDUP (
	SAMTOOLS_INDEX.out.bambai
    )
    ch_software_versions = ch_software_versions.mix(UMITOOLS_DEDUP.out.version.ifEmpty(null))

    ch_input_unencoded = Channel.from('G')
    UNENCODED_G (
	UMITOOLS_DEDUP.out.bam,
	ch_input_unencoded
    )
    ch_software_versions = ch_software_versions.mix(UNENCODED_G.out.version.ifEmpty(null))

    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }


    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowReaptec.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))


  
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
