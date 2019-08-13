params.outdir = 'results'

include 'index' params(params)
include 'quant' params(params)
include 'fastqc' params(params)
include 'multiqc' params(params)

workflow RNASEQ {
  get:
    transcriptome
    read_pairs_ch
    multiqc_config
 
  main: 
    INDEX(transcriptome)

    FASTQC(read_pairs_ch)

    QUANT(INDEX.out, read_pairs_ch)

    MULTIQC( 
        QUANT.out.mix(FASTQC.out).collect(),
        multiqc_config )

  emit: 
     QUANT.out
}