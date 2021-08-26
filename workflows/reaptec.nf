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
			  params.barcodes,
                          params.ref_chrom,
                          params.ref_pro1,
                          params.ref_enh,
			  params.enhancers_mask
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.barcodes) { ch_barcodes = file(params.barcodes) } else { exit 1, 'Input barcodes not specified!' }
if (params.ref_chrom) { ch_chrom = file(params.ref_chrom) } else { exit 1, 'Input reference chromosome sizes not specified!' }
if (params.ref_pro1) { ch_pro1 = file(params.ref_pro1) } else { exit 1, 'Input FANTOM5 promoters not specified!' }
if (params.ref_enh) { ch_enh = file(params.ref_enh) } else { exit 1, 'Input FANTOM5 enhancers not specified!' }
if (params.enhancers_mask) { ch_enhancers_mask = file(params.enhancers_mask) } else { exit 1, 'Input enhancers mask not specified!' }
if (params.enhancers_prefix) { ch_enhancers_prefix = params.enhancers_prefix } else { exit 1, 'Input prefix for the final call of the enhancers not specified!' }
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
def cat_fastq_options          = modules['cat_fastq']
if (!params.save_merged_fastq) { cat_fastq_options['publish_files'] = false }

params.star_index_options   = [:]

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC  }             from '../modules/nf-core/modules/fastqc/main'              addParams( options: modules['fastqc'] )
include { CAT_FASTQ}            from '../modules/nf-core/modules/cat/fastq/main'           addParams( options: cat_fastq_options)
include { MULTIQC }             from '../modules/nf-core/modules/multiqc/main'             addParams( options: multiqc_options   )
include { UMITOOLS_EXTRACT }    from '../modules/nf-core/modules/umitools/extract/main'    addParams( options: modules['umitools_extract'])
include { CUTADAPT }            from '../modules/nf-core/modules/cutadapt/main'            addParams( options: modules['cutadapt'] )
include { GET_WHITELIST }       from '../modules/local/get_whitelist'                      addParams( options: [:] )
include { STAR_GENOMEGENERATE } from '../modules/nf-core/modules/star/genomegenerate/main' addParams( options: modules['star_genomegenerate'])
include { STAR_ALIGN }          from '../modules/nf-core/modules/star/align/main'          addParams( options: modules['star_align'])
include { SAMTOOLS_INDEX }      from '../modules/nf-core/modules/samtools/index/main'      addParams( options: [:])
include { UMITOOLS_DEDUP }      from '../modules/nf-core/modules/umitools/dedup/main'      addParams( options: modules['umitools_dedup'])
include { UNENCODED_G }         from '../modules/local/unencoded_g'                        addParams( options: [:] )
include { BAM_TO_CTSS }         from '../modules/local/bam_to_ctss'                        addParams( options: modules['bam_to_ctss'] )
include { CALL_ENHANCERS }      from '../modules/local/call_enhancers'                     addParams( options: modules['call_enhancers'] )
include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GTF }      from '../modules/nf-core/modules/gunzip/main'              addParams( options: [:] )
include {
    UNTAR as UNTAR_STAR_INDEX } from '../modules/nf-core/modules/untar/main'               addParams( options: params.star_index_options   )

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
    .map {
            meta, fastq ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, fastq ] }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_fastq.multiple
    )
    ch_software_versions = ch_software_versions.mix(FASTQC.out.version.ifEmpty(null))


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
	ch_cat_fastq,
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
	.set { ch_fastq_r1 }

//3) Trim the TSO sequence (13 nt) with cutadapt from Read1 processed by umi-tools.
// #3) Trim the TSO sequence (13 nt) with cutadapt from Read1 processed by umi-tools.
// cutadapt -u 13 -o /local/home/ubuntu/DATA/ForIFOM/Process_Test_S1_L001_R1_001_trim13.fastq.gz /local/home/ubuntu/DATA/ForIFOM/Process_Test_S1_L001_R1_001.fastq.gz -j 20
    
    CUTADAPT (
	ch_fastq_r1
    )
    ch_software_versions = ch_software_versions.mix(CUTADAPT.out.version.ifEmpty(null))


// #4) Map read1s only using STAR
// STAR --runThreadN 12 --genomeDir /local/home/ubuntu/Ref/homSap_hg38_GENCODEv34/indexes_STAR-2.6.0c_homSap_GRCh38.p13_GENCODEv34_cmp_chr/ \
// --readFilesIn /local/home/ubuntu/DATA/ForIFOM/Process_Test_S1_L001_R1_001_trim13.fastq.gz \
// --readFilesCommand zcat \
// --outFilterMultimapNmax 1 --outTmpDir /local/home/ubuntu/DATA/ForIFOM/STAR_Tmp/ForIFOM_ \
    // --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /local/home/ubuntu/DATA/ForIFOM/STAR_results/ForIFOM_

    // ref: https://github.com/nf-core/rnaseq/blob/master/subworkflows/local/prepare_genome.nf#L55
    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (params.gtf.endsWith('.gz')) {
        ch_genome_gtf = GUNZIP_GTF ( params.gtf ).gunzip
    } else {
        ch_genome_gtf = file(params.gtf)
    }

    // Ref: https://github.com/nf-core/rnaseq/blob/master/subworkflows/local/prepare_genome.nf#L126
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index   = Channel.empty()
    ch_star_version = Channel.empty()
    if (params.star_index) {

        if (params.star_index.endsWith('.tar.gz')) {
            ch_star_index = UNTAR_STAR_INDEX ( params.star_index ).untar
        } else {
            ch_star_index = file(params.star_index)
        }
    } else {

	// ref: https://github.com/nf-core/rnaseq/blob/master/subworkflows/local/prepare_genome.nf#L46
	//
	// Uncompress genome fasta file if required
	//
	if (params.fasta.endsWith('.gz')) {
            ch_genome_fasta = GUNZIP_FASTA ( params.fasta ).gunzip
        } else {
            ch_genome_fasta = file(params.fasta)
        }

        ch_star_index   = STAR_GENOMEGENERATE ( ch_genome_fasta, ch_genome_gtf ).index
        ch_star_version = STAR_GENOMEGENERATE.out.version
    }

    ch_software_versions = ch_software_versions.mix(ch_star_version.ifEmpty(null))


    STAR_ALIGN (
	CUTADAPT.out.reads,
	ch_star_index,
	ch_genome_gtf
    )
    ch_software_versions = ch_software_versions.mix(STAR_ALIGN.out.version.ifEmpty(null))


//     #5) Deduplicate with umi-tools
// samtools index -@ 8 ForIFOM_Aligned.sortedByCoord.out.bam
// umi_tools dedup --per-cell -I ForIFOM_Aligned.sortedByCoord.out.bam --output-stats=deduplicated -S ForIFOM_deduplicated_Nopairedoption.bam

    SAMTOOLS_INDEX (
	STAR_ALIGN.out.bam
    )
    ch_software_versions = ch_software_versions.mix(SAMTOOLS_INDEX.out.version.ifEmpty(null))

    UMITOOLS_DEDUP (
	// combines the bam with its bai
	STAR_ALIGN.out.bam.join(SAMTOOLS_INDEX.out.bai, by: [0])
    )
    ch_software_versions = ch_software_versions.mix(UMITOOLS_DEDUP.out.version.ifEmpty(null))


//     #6) Extract reads that start with unencoded-G. (UnencodedG.sh)
// mkdir ForSfclG
// mv ForIFOM_deduplicated_Nopairedoption.bam ./ForSfclG
// cd ForSfclG
// /local/home/ubuntu/Scripts/UnencodedG.sh G

    ch_input_unencoded = Channel.from('G')
    UNENCODED_G (
	UMITOOLS_DEDUP.out.bam,
	ch_input_unencoded
    )
    ch_software_versions = ch_software_versions.mix(UNENCODED_G.out.version.ifEmpty(null))


    //TODO: implement the filtering on cell clusters providing the list of barcodes to select
    
// 8) Convert the unencodedG.bam file to the ctss file. (BAMtoCTSS.sh)
// (+ we can check how much reads overlap FANTOM5 promoters and FANTOM5 enhancers)

// ## bedGraphToBigWig (http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/)
// ## bigWigAverageOverBed (http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/)
// ## bigWigToBedGraph (http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/)
// ## bigWigMerge (http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/)
// ## samtools: 1.10 (I use this version)
// ## bedtools: v2.29.2 (I use this version)

// ex.)(please change the PATH to the working directory and STAR index)
// ## Please change the PATH to Ref_chrom,Ref_pro1 and Ref_enh in the BAMtoCTSS.sh.

// cd /local/home/ubuntu/DATA/ForIFOM/STAR_results/ForSfclG
// mkdir SfclG
// mv SoftclipG* ./SfclG
// cd SfclG
// /local/home/ubuntu/Scripts/BAMtoCTSS.sh


    BAM_TO_CTSS (
	ch_star_index,
	ch_pro1,
	ch_enh,
	UNENCODED_G.out.clipped
    )
    ch_software_versions = ch_software_versions.mix(BAM_TO_CTSS.out.version.ifEmpty(null))

// 9-2) Call all enhancers.

// cd /local/home/ubuntu/DATA/ForIFOM/STAR_results/ForSfclG/SfclG/SoftclipG_ForIFOM_deduplicated_Nopairedoption
// gunzip *ctss.bed.gz
// mkdir enhancer_300bpmask_RESULTs
// find `pwd` -name SoftclipG_ForIFOM_deduplicated_Nopairedoption_mq20.ctss.bed > SoftclipG_ForIFOM_deduplicated_Nopairedoption_mq20.ctss.txt

// cd /local/home/ubuntu/enhancers/scripts/
// CTSS='/local/home/ubuntu/DATA/ForIFOM/STAR_results/ForSfclG/SfclG/SoftclipG_ForIFOM_deduplicated_Nopairedoption'
// MASK='/local/home/ubuntu/Ref/gencode.v34.transcript.ptrcoding.300bpslop.mask.bed'
// ./fixed_bidir_enhancers_10bp -s ForIFOM_all_0.8_0_10bp_ -m ${MASK} -t ${CTSS}/enhancer_300bpmask_RESULTs/Tmp_ForIFOM_all_0.8_0_10bp -f ${CTSS}/*ctss.txt -o ${CTSS}/enhancer_300bpmask_RESULTs/Result_ForIFOM_all_0.8_0_10bp

// ## you can get ./enhancer_300bpmask_RESULTs/Result_ForIFOM_all_0.8_0_10bp/ForIFOM_all_0.8_0_10bp_enhancers.bed
// ## ForIFOM_all_0.8_0_10bp_enhancers.bed is enhancer bed file detected from the sample. (This time 2,058 enhancers.)

    CALL_ENHANCERS (
	BAM_TO_CTSS.out.ctss,
	ch_enhancers_mask,
	ch_enhancers_prefix
    )
    ch_software_versions = ch_software_versions.mix(CALL_ENHANCERS.out.version.ifEmpty(null))




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
