#!/usr/bin/env nextflow 

processDef index {
    tag "$transcriptome_file.simpleName"

    input:
    file transcriptome 

    output:
    file 'index' 

    script:
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}


processDef quant {
    tag "$pair_id"

    input:
    file index 
    set pair_id, file(reads) 

    output:
    file(pair_id) 

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 ${reads[0]} -2 ${reads[1]} -o $pair_id
    """
}

processDef fastqc {
    tag "FASTQC on $sample_id"
    publishDir params.outdir

    input:
    set sample_id, file(reads)

    output:
    file("fastqc_${sample_id}_logs") 

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}


processDef multiqc {
    publishDir params.outdir, mode:'copy'

    input:
    file('*') 
    file(config) 

    output:
    file('multiqc_report.html')

    script:
    """
    cp $config/* .
    echo "custom_logo: \$PWD/logo.png" >> multiqc_config.yaml
    multiqc .
    """
}
