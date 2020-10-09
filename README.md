# VizFaDa RNA-seq quant pipeline

This Nextflow pipeline allows fast transcript quantification of short-read RNA-seq data from [FAANG](https://data.faang.org).

It uses [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) for (optional) quality control,
 [Salmon](https://salmon.readthedocs.io/en/latest/salmon.html) for transcript quantification,
 and [tximport](http://bioconductor.org/packages/release/bioc/html/tximport.html) to output transcript and gene level quantification matrices.

It can pull fastq files from the FAANG data portal either:
* by FIRE access if working from the EMBL-EBI [Embassy Cloud](https://www.embassycloud.org/)
* using [Aspera downloads](https://www.ibm.com/products/aspera/downloads) from the [EMBL-EBI ENA](https://www.ebi.ac.uk/ena/browser/home)
* uswing wget from the [EMBL-EBI ENA](https://www.ebi.ac.uk/ena/browser/home)

## Example usage

```bash
nextflow run lauramble/rnaseq-vizfada \
    -r v2.0 \
    -resume \
    --species Gallus_gallus \
    --all true \
    --data ''
```

## Pipeline parameters

* `all`: boolean (`true`/`false`, default: `false`). Should all available FAANG RNA-seq matching `species` be analysed, or only the subset specified at `input`?

Example:

```bash
nextflow run lauramble/rnaseq-vizfada \
    -r v2.0 \
    --species Gallus_gallus \
    --all true \
    --data ''
```

* `input`: file path (default: `"$baseDir/data/test_input.txt"`). If `all` is `false`, the pipeline will analyse only the FAANG RNAseq specified in the file.
The `input` file should be a text file containing one ENA Run ID per line, from FAANG experiments, i.e. :
```
ERR1464185

```

* `data`: 

* `fastqc`

* `input`

* `keepReads`

* `species`

* `species_ensembl`

* `aspera`

* `asperaPath`

* `fire`

* Other parameters: `custom_config_version`, custom_config_base, max_cpus, max_memory, multiqc, outdir, salmon



## Requirements

* Unix-like operating system (Linux, macOS, etc)
* Java 8

## Quickstart

1. If you don't have it already install Docker in your computer. Read more [here](https://docs.docker.com/).

2. Install Nextflow (version 0.24.x or higher):

        curl -s https://get.nextflow.io | bash

3. Launch the pipeline execution:

        ./nextflow run nextflow-io/rnaseq-nf -with-docker

4. When the execution completes open in your browser the report generated at the following path:

        results/multiqc_report.html

You can see an example report at the following [link](http://multiqc.info/examples/rna-seq/multiqc_report.html).

Note: the very first time you execute it, it will take a few minutes to download the pipeline
from this GitHub repository and the the associated Docker images needed to execute the pipeline.  


## Cluster support

RNASeq-NF execution relies on [Nextflow](http://www.nextflow.io) framework which provides an
abstraction between the pipeline functional logic and the underlying processing system.

This allows the execution of the pipeline in a single computer or in a HPC cluster without modifying it.

Currently the following resource manager platforms are supported:

  + Univa Grid Engine (UGE)
  + Platform LSF
  + SLURM
  + PBS/Torque


By default the pipeline is parallelized by spawning multiple threads in the machine where the script is launched.

To submit the execution to a UGE cluster create a file named `nextflow.config` in the directory
where the pipeline is going to be executed with the following content:

    process {
      executor='uge'
      queue='<queue name>'
    }

To lean more about the avaible settings and the configuration file read the
Nextflow [documentation](http://www.nextflow.io/docs/latest/config.html).


## Components

RNASeq-NF uses the following software components and tools:

* [Salmon](https://combine-lab.github.io/salmon/) 0.8.2
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 0.11.5
* [Multiqc](https://multiqc.info) 1.0
