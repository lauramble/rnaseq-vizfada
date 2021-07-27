/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// Get list of allowed species
def valid_params = [
    species : ["bos_taurus", "capra_hircus", "equus_caballus", "gallus_gallus", "ovis_aries", "sus_scrofa"],
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRnaseqVizFaDa.initialise(params, log, valid_params)

SPECIES = "${params.species}".toLowerCase().replace(" ", "_")

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

// RNASEQ SUBWORKFLOWS

def prepareToolIndices  = 'salmon'
def biotype = 'gen_biotype'

def salmon_index_options   = modules['salmon_index']
def salmon_quant_options   = modules['salmon_quant']
def rsem_preparereference_options = modules['rsem_preparereference']


//include { INPUT_CHECK    } from '../subworkflows/local/input_check'    addParams( options: [:] )
/*
include { GET_FAANG } from '../modules/local/get_faang'
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome' addParams(outdir: "${params.index}/${SPECIES}",
                                                                                 rsem_index_options: rsem_preparereference_options,
                                                                                 salmon_index_options: salmon_index_options )
include { FETCH_AND_RNASEQ } from './fetch_and_rnaseq' addParams( input: params.input, index: "${params.index}/${SPECIES}" )
*/
include { FETCH_AND_RNASEQ } from './fetch_and_rnaseq' addParams(index: "${params.index}" ) //TODO: add species

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow RNASEQ_VIZFADA {
    /*
    if ( !file("${params.index}/${SPECIES}").exists() )
      PREPARE_GENOME( prepareToolIndices, biotype )
    */
    
    println "=== RNASEQ_VIZFADA ==="
    
    ch_ids = Channel.fromPath(params.ids)
              
    FETCH_AND_RNASEQ( ch_ids )

}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    //NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    NfcoreTemplate.summary(workflow, params, log)
    //WorkflowFetchngs.curateSamplesheetWarn(log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
