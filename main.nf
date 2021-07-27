/*
========================================================================================
    VALIDATE INPUTS
========================================================================================


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

if (params.index) {
  file(params.index)
}


========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================


// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

// RNASEQ SUBWORKFLOWS

def salmon_index_options   = modules['salmon_index']
def salmon_quant_options   = modules['salmon_quant']

//include { INPUT_CHECK    } from '../subworkflows/local/input_check'    addParams( options: [:] )
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome' addParams( genome_options: publish_genome_options, index_options: publish_index_options, gffread_options: gffread_options,  star_index_options: star_genomegenerate_options,  hisat2_index_options: hisat2_build_options, rsem_index_options: rsem_preparereference_options, salmon_index_options: salmon_index_options )
include { FETCH_AND_RNASEQ } from '../subworkflows/local/fetch_and_rnaseq' addParams( )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================


workflow RNASEQ_VIZFADA {

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
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================


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

nextflow.enable.dsl = 2

include { RNASEQ_VIZFADA } from './workflows/rnaseq_vizfada'

//
// WORKFLOW: Run main nf-core/rnaseq analysis pipeline
//
workflow VIZFADA {
    RNASEQ_VIZFADA()
}

workflow {
    VIZFADA ()
}
