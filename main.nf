/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'RNASEQ-NF'.
 *
 *   RNASEQ-NF is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   RNASEQ-NF is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with RNASEQ-NF.  If not, see <http://www.gnu.org/licenses/>.
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

params.reads = "$baseDir/data/ggal/*_{1,2}.fq"
params.transcriptome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
params.outdir = "results"
params.multiqc = "$baseDir/multiqc"

log.info """\
 R N A S E Q - N F   P I P E L I N E
 ===================================
 transcriptome: ${params.transcriptome}
 reads        : ${params.reads}
 outdir       : ${params.outdir}
 """

// note: import must be after params def  
importLibrary 'rnaseq.nf'


transcriptome_file = file(params.transcriptome)
multiqc_file = file(params.multiqc)


Channel
    .fromFilePairs( params.reads, checkIfExists: true )
    .into { read_pairs_ch; read_pairs2_ch }

index(transcriptome_file)
    .quant(read_pairs_ch)
    .mix(fastqc(read_pairs2_ch))
    .collect()
    .multiqc(multiqc_file)

workflow.onComplete {
	println ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
