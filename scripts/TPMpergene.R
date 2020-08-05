# packages="BiocManager"
# biocPackages=c("AnnotationHub", "tximport", "ensembldb")
# 
# install.packages(setdiff(packages, rownames(installed.packages())))
# 
# BiocManager::install(setdiff(biocPackages, rownames(installed.packages())), ask=F)

library(AnnotationHub)
library(tximport)
library(ensembldb)

args= commandArgs(trailingOnly = T)

if (length(args)!=4){
  stop("Rscript TPMpergene.R <quantDir> <species> <version> <release>", call. = F)
}

quantDir=args[1]
species=args[2]
version=args[3]
release=args[4]

# Query AnnotationHub for reference genome in Ensembl database
a=AnnotationHub::AnnotationHub()

q=AnnotationHub::query(a, c("EnsDb", species, version, release))
# r=stringr::str_extract(q$tags, "[0-9]+")
# r=max(as.numeric(r))

# q=AnnotationHub::query(a, c("EnsDb", species, version, r))
galgalEnsDb=q[[1]]


# Build tx2gene file for tximport
galgalTx=ensembldb::transcripts(galgalEnsDb, return.type="DataFrame")
galgalTx2gene = galgalTx[,c("tx_id_version", "gene_id")] # /!\ Use tx_id_version

# Get run names and quant files
files=list.files(path=quantDir, pattern="quant.sf", full.names=T, recursive=T)
names=list.dirs(path=quantDir, full.names=F, recursive=F)

# Tximport
galgalTximport=tximport::tximport(files, type="salmon", tx2gene = galgalTx2gene)

# Get and write feature table
galgalTPM=galgalTximport$abundance
colnames(galgalTPM)=names
write.table(galgalTPM, "GeneMatrixTPM.tsv", sep="\t", quote = F, col.names = NA)

# gal_abundance=read.table("abundance.csv", sep=";", header=T)
# row.names(gal_abundance)=gal_abundance$X
# galgalTPM=gal_abundance[,-1]

# Create heatmap
# galgalCor=cor(log10(galgalTPM+1))
# heatmap(galgalCor, symm=T)
