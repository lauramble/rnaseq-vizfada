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
 * Proof of concept of a RNAseq pipeline implemented with Nextflow
 *
 * Authors:
 * - Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * - Emilio Palumbo <emiliopalumbo@gmail.com>
 * - Evan Floden <evanfloden@gmail.com>
 */


/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.species_ensembl = "$baseDir/data/species_ensembl.txt"
params.outdir = "results"
params.multiqc = "$baseDir/multiqc"
params.fastqc = false
params.salmon = ""
params.species = "Gallus gallus"
params.input = "$baseDir/data/test_input.txt"
params.cpus = 3
params.data = "/data"

species=params.species.replaceAll(/ /, "_")
index=file("${params.data}/${species}/index", type:"dir").exists()

log.info """\
 R N A S E Q - N F   P I P E L I N E
 ===================================
 species      : ${params.species}
 baseDir      : ${baseDir}
 index        : ${index}
 outdir       : ${params.outdir}
 cpus         : ${params.cpus}
 fastqc       : ${params.fastqc}
 input        : ${params.input}
 """

if (!index) {
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
        publishDir "$params.data/$species", mode:'copy'
        publishDir params.outdir
        cpus params.cpus

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
    index_ch = Channel.fromPath( "${params.data}/${species}/index" )
}
 
 
ch_input=Channel.fromList(file(params.input).readLines())

process dlFromFaang {
    tag "$accession"
    // maxForks 1

    input:
    each accession from ch_input

    output:
    tuple val(accession), file("${accession}_1.fastq.gz"), file("${accession}_2.fastq.gz") into read_pairs_ch, read_pairs2_ch
    
    shell:
    if (params.fire) {
        """
        for accession in ${accession}_1 ${accession}_2
        do 
          url=\$(wget "http://data.faang.org/api/file/\$accession" -q -O - | grep -Po "/fastq/.*\\.fastq\\.gz")
          url=https://hh.fire.sdo.ebi.ac.uk/fire/public/era\$url
          wget \$url
        done
        """
    } else {
        """
        for accession in ${accession}_1 ${accession}_2
        do 
          url=\$(wget "http://data.faang.org/api/file/\$accession" -q -O - | grep -Po "ftp.*\\.fastq\\.gz")
          wget \$url
        done
        """
    }
}

/*
if ( params.index != "" ) {
    index_ch = Channel.fromPath( params.index )
} else {
    process index {
        tag "$transcriptome.simpleName"
        publishDir params.outdir, mode:'copy'
        cpus params.cpus

        input:
        path transcriptome from params.transcriptome

        output:
        path 'index' into index_ch

        script:
        """
        # code some changes in your script and save them
        salmon index --threads $task.cpus -t $transcriptome -i index
        """
    }
}
*/

process quant {
    tag "$pair_id"
    publishDir "${params.outdir}/quant", mode:'copy'
    cpus params.cpus

    input:
    path index from index_ch
    tuple val(pair_id), path(reads_1), path(reads_2) from read_pairs_ch

    output:
    path(pair_id) into (quant_ch, quant2_ch)

    script:
    """
    salmon quant --threads $task.cpus --libType=U -i $index -1 $reads_1 -2 $reads_2 -o $pair_id $params.salmon
    rm -f $reads_1 $reads_2
    """
}

if ( params.fastqc ) {

    process fastqc {
        tag "FASTQC on $sample_id"
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(sample_id), path(reads) from read_pairs2_ch

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
    container 'jmeigs1/rscript'
    publishDir params.outdir, mode:'copy'
    
    input:
    path "quant" from quant2_ch.collect()
    
    output:
    file abundance
    
    shell:
    """
    Rscript $baseDir/scripts/TPMpergene.R $outdir/quant
    """
}

workflow.onComplete {
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}

