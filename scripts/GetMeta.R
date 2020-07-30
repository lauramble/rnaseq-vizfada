library(jsonlite)
library(dplyr)

args= commandArgs(trailingOnly = T)

if (length(args)!=3){
  stop("Rscript GetMeta.R <path-to-specimens> <path-to-experiments> <path-to-species>", call. = F)
}

specimens=args[1]
experiments=args[2]
species=args[3]

#' Parse FAANG Metadata (obtained with FAANG's ES API) and add relevant specimen data.
#' 
#' @param metaFAANG A JSON file returned by FAANG's search API. Must be the result of a search on the 'files' index.
#' @param specimens A JSON file containing specimens metadata
#'
#' @return A data.frame
parseMetaFAANG = function(metaFAANG, specimens="../data/MetaFAANG/specimens_faang.json"){
  metaDf=jsonlite::read_json(metaFAANG, simplifyVector=T)
  metaDf=jsonlite::flatten(metaDf$hits$hits)
  
  specimens=jsonlite::read_json(specimens, simplifyVector = T)
  specimensDf=jsonlite::flatten(specimens$hits$hits)
  
  allDf=merge(x=metaDf, y=specimensDf, by.x="_source.specimen", by.y="_id", all.x=T)
  return(allDf)
}

#' Filter parsed FAANG metadata (with \link{parseMetaFAANG}) with a JSON file of experiments returned by FAANG's search API.
#' 
#' @param metaDF 
#' @param experiments 
#'
#' @return A data.frame
filterExperiments = function(metaDF, experiments="../data/MetaFAANG/experiments_faang.json") {
  experiments=jsonlite::read_json(experiments, simplifyVector=T)
  experimentsDf=jsonlite::flatten(experiments$hits$hits)
  
  filtMetaDf=metaDF[metaDF$`_source.experiment.accession` %in% experimentsDf$`_id`,]
  return(filtMetaDf)
}

### Parse FAANG metadata and filter only RNA-Seq experiments

# Use files obtained with the 'extraction_faang.sh' script
galgal=parseMetaFAANG(species, specimens = specimens)
galgalRNA=filterExperiments(galgal, experiments = experiments)
galgalExport=apply(galgalRNA,2,as.character) # For csv export

cat(unique(galgalRNA$`_source.run.accession`), sep="\n", file="input.txt")
write.table(galgalExport, file="metadata.tsv", sep="\t", row.names=F)

