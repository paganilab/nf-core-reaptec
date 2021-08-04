#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/reaptec
========================================================================================
    Github : https://github.com/nf-core/reaptec
    Website: https://nf-co.re/reaptec
    Slack  : https://nfcore.slack.com/channels/reaptec
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.gtf = WorkflowMain.getGenomeAttribute(params, 'gtf')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { REAPTEC } from './workflows/reaptec'

//
// WORKFLOW: Run main nf-core/reaptec analysis pipeline
//
workflow NFCORE_REAPTEC {
    REAPTEC ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_REAPTEC ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
