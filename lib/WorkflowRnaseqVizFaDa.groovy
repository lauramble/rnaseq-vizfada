//
// This file holds several functions specific to the workflow/rnaseq.nf in the nf-core/rnaseq pipeline
//

class WorkflowRnaseqVizFaDa {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log, valid_params) {

        if (!params.species) {
            log.error "No species provided. Check the file data/species_ensembl.tsv to get available species."
            System.exit(1)
        }

        if (!params.input && !params.all) {
            log.error "No parameter --input <FILE.txt> or --all. Provide an input file (check documentation) or use the parameter --all to process all the available RNASeq data from FAANG."
            System.exit(1)
        }
        
        if (!valid_params['species'].contains(params.species)){
          log.error "Invalid species: '${params.species}'. Available species: ${valid_params['species'].join(', ')}."
        }

        if (params.all) {
            if (params.input) {
                allInputWarn(log)
            }
        }
    }
    
    //
    // Print a warning if both all and input have been provided
    //
    private static void allInputWarn(log) {
        log.warn "=============================================================================\n" +
            "  Both '--all' and '--input' parameters have been provided.\n" +
            "  Using input file as priority.\n" +
            "==================================================================================="
    }
}
