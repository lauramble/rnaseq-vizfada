/*
 * Copyright (c) 2013-2019, Centre for Genomic Regulation (CRG).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
 *
 */


/*
 * RNASeq pipeline for VizFaDa
 *
 * Authors:
 * - Laura Morel <laura.morel@inrae.fr>
 * - Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * - Emilio Palumbo <emiliopalumbo@gmail.com>
 * - Evan Floden <evanfloden@gmail.com>
 */


/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */


species=params.species.replaceAll(/ /, "_")
index=file("${params.data}/${species}/index", type:"dir")

log.info """\
 R N A S E Q - N F   P I P E L I N E
 ===================================
 species      : ${params.species}
 baseDir      : ${baseDir}
 index        : ${index}
 outdir       : ${params.outdir}
 fastqc       : ${params.fastqc}
 input        : ${params.input}
 """

if (params.all) {
    process getMetaAndInput {
        tag "$species"
        label "R"
        
        container 'lauramble/r-vizfada:latest'
        publishDir "$params.outdir", pattern: 'metadata.tsv', mode: 'copy'
        
        input:
        val species from Channel.from(params.species)
        file 'extraction_faang.sh' from Channel.fromPath("$baseDir/scripts/extraction_faang.sh")
        file 'GetMeta.R' from Channel.fromPath("$baseDir/scripts/GetMeta.R")
        
        output:
        file 'input.txt' into input
        file 'metadata.tsv'
        
        shell:
        """
        bash extraction_faang.sh '$species' &> temp.txt
        Rscript GetMeta.R specimens.json experiments.json species.json
        """      
    }
    ch_input=input.map{ it -> it.readLines() }.flatten()
} else {
    ch_input=Channel.fromList(file(params.input).readLines())
}

if (!index.exists()) {
    process getcDNA {
        tag "$species"
        
        input:
        val species
        
        output:
        file "*.fa.gz" into ch_transcriptome
        
        shell:
        """
        version=\$( grep $species $params.species_ensembl | awk '{print \$2}' )
        url=\$( awk 'NR==1{print \$2}' $params.species_ensembl | sed "s/SPECIES/\\L$species/1; s/SPECIES.VERSION/$species.\$version/1" )
        wget \$url
        """
    }
    
    process index {
        tag "$transcriptome.simpleName"
        label "salmon"
        
        publishDir "$params.data/$species", mode:'copy'
        publishDir params.outdir, mode:'copy'

        input:
        path transcriptome from ch_transcriptome

        output:
        path "index" into index_ch

        script:
        """
        # code some changes in your script and save them
        salmon index --threads $task.cpus -t $transcriptome -i index
        """
    }
} else {
    index_ch=Channel.fromPath(index)
}
 

if (params.fire){
    baseURL='https://hh.fire.sdo.ebi.ac.uk/fire/public/era'
    regex='"/fastq/.*\\.fastq\\.gz"'
} else if (params.aspera && params.asperaPath) {
    baseURL='-QT -P33001 -i /tools/aspera/etc/asperaweb_id_dsa.openssh era-fasp@fasp'
    regex='"\\.sra.*\\.fastq\\.gz"'
} else {
    baseURL=''
    regex='"ftp.*\\.fastq\\.gz"'
}

process dlFromFaangAndQuant {
    tag "$accession"
    label "salmon"
    label "canIgnore"
    
    if (params.keepReads) {publishDir "${params.outdir}/reads", pattern: "*.fastq.gz", mode: 'copy'}
    publishDir "${params.outdir}/quant", mode:'copy', pattern: "${accession}"
    
    input:
    each accession from ch_input
    path index from index_ch
    val asperaPath from Channel.from(params.asperaPath)

    output:
    path "${accession}" into quant_ch, quant2_ch
    
    shell:
    '''
    #!/bin/bash
    checkpaired=$(wget "http://data.faang.org/api/file/_search/?size=25000" --post-data '{"query": { "wildcard": {"name": "!{accession}*"}}}' -q -O - | grep -Po "!{accession}(_[12])+")
    
    if (( $(echo $checkpaired | wc -w) != 0 ))
    then
        files="!{accession}_1 !{accession}_2"
    else
        files=!{accession}
    fi
    echo $checkpaired
    echo $files
    
    md5=""
    for file in $files
    do 
      url=$(wget "http://data.faang.org/api/file/$file" -q -O - | grep -Po !{regex})
      url=!{baseURL}$url
      checksum=$(wget http://data.faang.org/api/file/$file -q -O - | grep '"checksum": ".*?",' -Po | cut -d'"' -f4)
      if [[ !{params.aspera} == "true" ]]
      then
        !{asperaPath} !{baseURL}$url
      else
        while [[ $md5 != $checksum ]]
        do
            wget $url -O $file".fastq.gz"
            md5=$(md5sum $file".fastq.gz" | cut -d" " -f1)
        done
      fi
    done
    
    if (( $(echo $files | wc -w) == 2))
    then
        reads_1=!{accession}_1.fastq.gz
        reads_2=!{accession}_2.fastq.gz
        salmon quant --threads !{task.cpus} -l A -i !{index} -1 $reads_1 -2 $reads_2 -o !{accession} !{params.salmon}
    else
        salmon quant --threads !{task.cpus} -l A -i !{index} -r !{accession}.fastq.gz -o !{accession} !{params.salmon}
    fi
    
    rm -rf *.fastq.gz
    '''
}

if ( params.fastqc ) {

    process fastqc {
        tag "$sample_id"
        publishDir "${params.outdir}/FastQC", mode:'copy'

        input:
        tuple val(sample_id), path(reads) from reads_for_fastqc_ch
        output:
        path "fastqc_${sample_id}_logs" into fastqc_ch

        script:
        """
        mkdir fastqc_${sample_id}_logs
        fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
        """
    }

} else {
    fastqc_ch = Channel.empty()
}

process multiqc {
    publishDir params.outdir, mode:'copy'
    
    label "canIgnore"
    
    input:
    path 'data*/*' from quant_ch.mix(fastqc_ch).collect()
    path config from params.multiqc

    output:
    path 'multiqc_report.html'

    script:
    """
    cp $config/* .
    echo "custom_logo: \$PWD/logo.png" >> multiqc_config.yaml
    multiqc -v .
    """
}

process tximport {
    label "R"
    label "canIgnore"
    
    container 'lauramble/r-vizfada:latest'
    publishDir params.outdir, mode:'copy'
    
    input:
    path "quant" from quant2_ch.collect()
    
    output:
    file "abundance.csv"
    
    script:
    """
    version=\$( grep $species $params.species_ensembl | awk '{print \$2}' )
    Rscript $baseDir/scripts/TPMpergene.R . "${params.species}" \$version
    """
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}

