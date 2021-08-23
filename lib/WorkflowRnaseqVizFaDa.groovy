//
// This file holds several functions specific to the workflow/rnaseq.nf in the nf-core/rnaseq pipeline
//

class WorkflowRnaseqVizFaDa {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log, valid_params) {

        if (!params.species && !params.index) {
            log.error "No species or index provided."
            System.exit(1)
        }

        if (!params.ids) {
            log.warn "No parameter --ids <FILE.txt>. Processing all available data for species ${params.species}."
        }
        
        if (!valid_params['species'].contains("${params.species}".toLowerCase().replace(" ", "_"))) {
          log.warn "Invalid species: '${params.species}'. Available species: ${valid_params['species'].join(', ')}."
        }
        
    }
}
