ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def cat_fastq_options          = modules['cat_fastq']


include { CAT_FASTQ             } from '../modules/nf-core/modules/cat/fastq/main'             addParams( options: cat_fastq_options                            )
include { FASTQC_UMITOOLS_TRIMGALORE } from '../subworkflows/nf-core/fastqc_umitools_trimgalore' addParams( fastqc_options: modules['fastqc'] )
include { QUANTIFY_SALMON } from '../subworkflows/local/quantify_salmon'    addParams( tximport_options: modules['salmon_tximport'], salmon_quant_options: salmon_quant_options, merge_counts_options: modules['salmon_merge_counts'] )


include { SRA_IDS_TO_RUNINFO      } from '../modules/local/sra_ids_to_runinfo'      addParams( options: modules['sra_ids_to_runinfo']      )
include { SRA_RUNINFO_TO_FTP      } from '../modules/local/sra_runinfo_to_ftp'      addParams( options: modules['sra_runinfo_to_ftp']      )
include { SRA_FASTQ_FTP           } from '../modules/local/sra_fastq_ftp'           addParams( options: modules['sra_fastq_ftp']           )
include { SRA_TO_SAMPLESHEET      } from '../modules/local/sra_to_samplesheet'      addParams( options: modules['sra_to_samplesheet'], results_dir: modules['sra_fastq_ftp'].publish_dir )
include { SRA_MERGE_SAMPLESHEET   } from '../modules/local/sra_merge_samplesheet'   addParams( options: modules['sra_merge_samplesheet']   )
include { MULTIQC_MAPPINGS_CONFIG } from '../modules/local/multiqc_mappings_config' addParams( options: modules['multiqc_mappings_config'] )
include { GET_SOFTWARE_VERSIONS   } from '../modules/local/get_software_versions'   addParams( options: [publish_files : ['tsv':'']]       )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow FETCHNGS {
    take:
    input //path to input file
    index //path to salmon index

    ch_software_versions = Channel.empty()

    //
    // MODULE: Get SRA run information for public database ids
    //
    SRA_IDS_TO_RUNINFO (
        ch_ids,
        params.ena_metadata_fields ?: ''
    )

    //
    // MODULE: Parse SRA run information, create file containing FTP links and read into workflow as [ meta, [reads] ]
    //
    SRA_RUNINFO_TO_FTP (
        SRA_IDS_TO_RUNINFO.out.tsv
    )

    SRA_RUNINFO_TO_FTP
        .out
        .tsv
        .splitCsv(header:true, sep:'\t')
        .map {
            meta ->
                meta.single_end = meta.single_end.toBoolean()
                [ meta, [ meta.fastq_1, meta.fastq_2 ] ]
        }
        .unique()
        .set { ch_sra_reads }
    ch_software_versions = ch_software_versions.mix(SRA_RUNINFO_TO_FTP.out.version.first().ifEmpty(null))

    if (!params.skip_fastq_download) {
        //
        // MODULE: If FTP link is provided in run information then download FastQ directly via FTP and validate with md5sums
        //
        SRA_FASTQ_FTP (
            ch_sra_reads.map { meta, reads -> if (meta.fastq_1)  [ meta, reads ] }
        )

        //
        // MODULE: Stage FastQ files downloaded by SRA together and auto-create a samplesheet
        //
        SRA_TO_SAMPLESHEET (
            SRA_FASTQ_FTP.out.fastq,
            params.nf_core_pipeline ?: '',
            params.sample_mapping_fields
        )

        //
        // MODULE: Create a merged samplesheet across all samples for the pipeline
        //
        SRA_MERGE_SAMPLESHEET (
            SRA_TO_SAMPLESHEET.out.samplesheet.collect{it[1]},
            SRA_TO_SAMPLESHEET.out.mappings.collect{it[1]}
        )

        //
        // MODULE: Create a MutiQC config file with sample name mappings
        //
        if (params.sample_mapping_fields) {
            MULTIQC_MAPPINGS_CONFIG (
                SRA_MERGE_SAMPLESHEET.out.mappings
            )
        }

        //
        // If ids don't have a direct FTP download link write them to file for download outside of the pipeline
        //
        def no_ids_file = ["${params.outdir}", "${modules['sra_fastq_ftp'].publish_dir}", "IDS_NOT_DOWNLOADED.txt" ].join(File.separator)
        ch_sra_reads
            .map { meta, reads -> if (!meta.fastq_1) "${meta.id.split('_')[0..-2].join('_')}" }
            .unique()
            .collectFile(name: no_ids_file, sort: true, newLine: true)
    }
    
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    .map {
        meta, fastq ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, fastq ] }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters
    //
    FASTQC_UMITOOLS_TRIMGALORE (
        ch_cat_fastq,
        params.skip_fastqc || params.skip_qc,
        params.with_umi,
        params.skip_trimming
    )
    ch_software_versions = ch_software_versions.mix(FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(FASTQC_UMITOOLS_TRIMGALORE.out.umitools_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(FASTQC_UMITOOLS_TRIMGALORE.out.trimgalore_version.first().ifEmpty(null))

    //
    // SUBWORKFLOW: Pseudo-alignment and quantification with Salmon
    //
    ch_salmon_multiqc                   = Channel.empty()
    ch_pseudoaligner_pca_multiqc        = Channel.empty()
    ch_pseudoaligner_clustering_multiqc = Channel.empty()
    QUANTIFY_SALMON (
        ch_trimmed_reads,
        params.salmon_index,
        ch_dummy_file,
        params.gtf,
        false,
        params.salmon_quant_libtype ?: ''
    )
    ch_salmon_multiqc = QUANTIFY_SALMON.out.results
    if (params.skip_alignment && params.aligner != 'star_salmon') {
        ch_software_versions = ch_software_versions.mix(QUANTIFY_SALMON.out.salmon_version.first().ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(QUANTIFY_SALMON.out.tximeta_version.ifEmpty(null))
        ch_software_versions = ch_software_versions.mix(QUANTIFY_SALMON.out.summarizedexperiment_version.ifEmpty(null))
    }

    if (!params.skip_qc & !params.skip_deseq2_qc) {
        DESEQ2_QC_SALMON (
            QUANTIFY_SALMON.out.counts_gene_length_scaled,
            ch_pca_header_multiqc,
            ch_clustering_header_multiqc
        )
        ch_pseudoaligner_pca_multiqc        = DESEQ2_QC_SALMON.out.pca_multiqc
        ch_pseudoaligner_clustering_multiqc = DESEQ2_QC_SALMON.out.dists_multiqc
        if (params.skip_alignment) {
            ch_software_versions = ch_software_versions.mix(DESEQ2_QC_SALMON.out.version.ifEmpty(null))
        }
    }

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
