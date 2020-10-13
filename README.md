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

* `species`: string (default: `Gallus gallus`).  Which FAANG species is analysed? should be one of the species listed
in [data/species_ensembl.txt](/data/species_ensembl.txt).

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
ERR1464184
ERR1464185
ERR1464186

```

* `fastqc`: boolean (`true`/`false`, default: `false`). Whether fastQC should be performed on the fastq files or not.
Running fastqc will make the pipeline slower.

* `keepReads`: boolean (`true`/`false`, default: `true`). Should the fastq files be kept or deleted after transcript quantification steps? Beware of potentially high disk usage if using `--all true --keepReads true`.

* `species_ensembl`: string (default to $baseDir/data/species_ensembl.txt).

* `aspera`: boolean (`true`/`false`, default: `false`). Should aspera download be used instead of wget? Aspera dowload is usually much faster. If `true`, `asperaPath` need tp be defined.

* `asperaPath`: path to the apsera binary. Example: `--asperaPath /tools/aspera/bin/ascp`

* `fire`: boolean (`true`/`false`, default: `false`). Should direct fire access be used instead of public ENA endpoints?
You will need to have fire access to the ENA data, most likely though an [EMBL-EBI Embassy Cloud instance](https://www.embassycloud.org/).

* `data`: name of the folder containing the fastq files. Default to `data`.

* `outdir`: path of the output directory. Default to `results`.

* Other parameters: `custom_config_version`, `custom_config_base`, `max_cpus`, `max_memory`, `multiqc`, `salmon`: see [nextflow.config](nextflow.config).


## Components

RNASeq-NF uses the following software components and tools:

* [Salmon](https://combine-lab.github.io/salmon/) 0.8.2
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 0.11.5
* [Multiqc](https://multiqc.info) 1.0
* [R](https://cran.r-project.org/)
* [AnnotationHub](http://bioconductor.org/packages/release/bioc/html/AnnotationHub.html)
* [tximport](http://bioconductor.org/packages/release/bioc/html/tximport.html)
* [ensembldb](http://bioconductor.org/packages/release/bioc/html/ensembldb.html)
