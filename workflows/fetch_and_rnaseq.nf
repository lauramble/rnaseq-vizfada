/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

//params.index = './index'
//def FASTA = file("${params.index}/*.fa")[0]
//def GTF = file("${params.index}/*.gtf")[0]

def GTF = file(params.gtf).baseName
def FASTA = file(params.fasta).baseName

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include { FETCHNGS } from '../subworkflows/local/fetchngs' addParams(nf_core_pipeline: 'rnaseq', outdir: "${params.outdir}/fetchngs")
include { RNASEQ } from '../subworkflows/local/rnaseq' addParams( skip_trimming: true,
                                                                  skip_alignment:true,
                                                                  pseudo_aligner: 'salmon',
                                                                  //salmon_index: "${params.index}/index/salmon", //TODO: change from index/salmon to salmon
                                                                  //fasta: "${params.index}/$FASTA",
                                                                  //gtf: "${params.index}/$GTF",
                                                                  //rsem_index: "${params.index}/rsem"
                                                                  )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow FETCH_AND_RNASEQ {
    take:
    pathToIDs //path to id list
    gtf
    salmon_index
    
    main:
    
    FETCHNGS( pathToIDs )
    
    RNASEQ( FETCHNGS.out.samplesheet, gtf, salmon_index )

}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/
/*
workflow.onComplete {
    NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    NfcoreTemplate.summary(workflow, params, log)
    WorkflowFetchngs.curateSamplesheetWarn(log)
}
*/
/*
========================================================================================
    THE END
========================================================================================
*/
