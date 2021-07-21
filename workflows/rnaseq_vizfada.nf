/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// Get list of allowed species
def allowed_species = []

file('species_ensembl.txt').withReader {
    String line
    while( line = it.readLine() ){
        allowed_species += line.split('\t')[0]
    }
}
allowed_species.remove(0)
allowed_species.remove(0)

def valid_params = [
    species : allowed_species,
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowRnaseqVizFaDa.initialise(params, log, valid_params)


/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

// RNASEQ SUBWORKFLOWS

def salmon_index_options   = modules['salmon_index']
def salmon_quant_options   = modules['salmon_quant']
def rsem_preparereference_options = modules['rsem_preparereference']


//include { INPUT_CHECK    } from '../subworkflows/local/input_check'    addParams( options: [:] )
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome' addParams(rsem_index_options: rsem_preparereference_options, salmon_index_options: salmon_index_options )
include { FETCH_AND_RNASEQ } from '../subworkflows/local/fetch_and_rnaseq' addParams( input: params.input, index: params.index )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow RNASEQ_VIZFADA {
    
    FETCH_AND_RNASEQ()
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
