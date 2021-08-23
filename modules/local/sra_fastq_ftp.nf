// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SRA_FASTQ_FTP {
    tag "$experiment"
    label 'process_medium'
    label 'error_retry'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:["id":experiment], publish_by_meta:['id']) }

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    path runinfo

    output:
    path "${experiment}_samplesheet.csv", emit: samplesheet
    //tuple val(meta), path("*md5")     , emit: md5

    shell:
    experiment = "${runinfo.getSimpleName()}"
    '''
    single_end=$(awk -F "\\t" '$2 !~/.*_1\\.fastq\\.gz$/ {if(NR>1) print $2}' !{runinfo})
    paired_end=$(awk -F "\\t" '$2 ~/.*_1\\.fastq\\.gz$/ {if(NR>1)print $2}' !{runinfo})

    cp !{runinfo} samplesheet.tsv

    for ftp in $paired_end
    do
    id=$(echo $ftp | sed 's/.*\\/\\(.*\\)\\/.*\\.fastq\\.gz/\\1/')
    ftp2=$(echo $ftp | sed s/_1/_2/)
    wget "ftp://$ftp" -O $id"_1.fastq.gz"
    wget "ftp://$(echo $ftp | sed s/_1/_2/)" -O $id"_2.fastq.gz"
    sed -i 0,/_T[0-9]*/{s/_T[0-9]*/_$id/} samplesheet.tsv
    sed -i "s~$ftp~$(readlink -f $id'_1.fastq.gz')~1" samplesheet.tsv
    sed -i "s~$ftp2~$(readlink -f $id'_2.fastq.gz')~1" samplesheet.tsv
    done

    for ftp in $single_end
    do
    id=$(echo $ftp | sed 's/.*\\/\\(.*\\)\\/.*\\.fastq\\.gz/\\1/')
    wget "ftp://$ftp" -O $id".fastq.gz"
    sed -i 0,/_T[0-9]*/{s/_T[0-9]*/_$id/} samplesheet.tsv
    sed -i "s~$ftp~$(readlink -f $id'.fastq.gz')~1" samplesheet.tsv
    done
    
    echo "samplesheet generated"
    cat samplesheet.tsv
    
    awk -F"\\t" 'BEGIN {OFS=","} {if(NR==1) print "sample",$2,$3,"strandedness",$1;else print "!{experiment}",$2,$3,"unstranded",$1}' samplesheet.tsv > !{experiment}_samplesheet.csv
    
    echo "final samplesheet"
    cat !{experiment}_samplesheet.csv
    '''

    /*
    script:
    if (meta.single_end) {
        """
        bash -c 'wget $options.args -L ftp://${fastq[0]} -O ${meta.id}.fastq.gz';

        echo "${meta.md5_1} ${meta.id}.fastq.gz" > ${meta.id}.fastq.gz.md5
        md5sum -c ${meta.id}.fastq.gz.md5
        """
    } else {
        """
        bash -c 'wget $options.args -L ftp://${fastq[0]} -O ${meta.id}_1.fastq.gz';

        echo "${meta.md5_1} ${meta.id}_1.fastq.gz" > ${meta.id}_1.fastq.gz.md5
        md5sum -c ${meta.id}_1.fastq.gz.md5

        bash -c 'wget $options.args -L ftp://${fastq[1]} -O ${meta.id}_2.fastq.gz';

        echo "${meta.md5_2} ${meta.id}_2.fastq.gz" > ${meta.id}_2.fastq.gz.md5
        md5sum -c ${meta.id}_2.fastq.gz.md5
        """
    }
    */
}
