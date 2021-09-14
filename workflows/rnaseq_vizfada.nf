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

params.fspecies = "${params.species}".toLowerCase().replace(" ", "_")
params.full_index = "${params.index}/${params.fspecies}"

def GTF = file(params.gtf).baseName
def FASTA = file(params.fasta).baseName

if (file(params.full_index).exists()){
  params.salmon_index = file("${params.full_index}/index/salmon").exists() ? file("${params.full_index}/index/salmon") : null
  params.rsem_index = file("${params.full_index}/rsem").exists() ? file("${params.full_index}/rsem") : null
  params.fasta = file("${params.full_index}/${FASTA}") ? file("${params.full_index}/${FASTA}") : params.fasta
  params.gtf = file("${params.full_index}/${GTF}") ? file("${params.full_index}/${GTF}") : params.gtf
}

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

salmon_index_options   = modules['salmon_index']
def salmon_quant_options   = modules['salmon_quant']
def rsem_preparereference_options = modules['rsem_preparereference']


//include { INPUT_CHECK    } from '../subworkflows/local/input_check'    addParams( options: [:] )

include { GET_FAANG } from '../modules/local/get_faang' addParams(ids: params.ids)
//include { FLOW_MANAGER } from '../modules/local/flow_manager'
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome' addParams(outdir: "${params.full_index}",
                                                                                 rsem_index_options: rsem_preparereference_options,
                                                                                 salmon_index_options: salmon_index_options )
//include { FETCH_AND_RNASEQ } from './fetch_and_rnaseq' addParams( input: params.input, index: "${params.index}/${params.fspecies}" )

include { FETCH_AND_RNASEQ } from './fetch_and_rnaseq' addParams(index: "${params.full_index}" ) //TODO: add species

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow RNASEQ_VIZFADA {

    PREPARE_GENOME( prepareToolIndices, biotype )
    
    GET_FAANG ()
    
    //ch_ids = Channel.fromPath(params.ids)
              
    FETCH_AND_RNASEQ( GET_FAANG.out.ids.flatten(),
                      PREPARE_GENOME.out.gtf,
                      PREPARE_GENOME.out.salmon_index)

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
