/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

params.index = './index'

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include { FETCHNGS } from '../subworkflows/local/fetchngs' addParams( input: params.input, nf_core_pipeline: 'rnasesq', outdir: './fetchngs')
include { RNASEQ } from '../subworkflows/local/rnaseq' addParams( skip_trimming: true,
                                                                  skip_alignment:true,
                                                                  pseudo_aligner: 'salmon',
                                                                  salmon_index: "${params.index}/index/salmon",
                                                                  rsem_index: "${params.index}/rsem")

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow FETCH_AND_RNASEQ {
  
    FETCHNGS()
    
    RNASEQ( FETCHNGS.out.samplesheet )
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    NfcoreTemplate.summary(workflow, params, log)
    WorkflowFetchngs.curateSamplesheetWarn(log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
