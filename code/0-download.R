#make sure to install these packages 
install.packages("BiocManager")
BiocManager::install(c( "GEOquery","clinfun", "GGally", "factoextra"))
# load libraries
library(tidyverse)
library(GEOquery)
library(canvasXpress)
library(clinfun)


  

# download gene expression and related info datasets from NIH GEO
data <- getGEO(GEO = "GSE112676")

# select related clinical data from NIH GEO
clindata <- data[["GSE112676_series_matrix.txt.gz"]]@phenoData@data

#select raw gene expression data from NIH GEO
gene_counts<- data[["GSE112676_series_matrix.txt.gz"]]@assayData[["exprs"]]

#keep rownames for gene data
gene_counts2=gene_counts%>%
  as_tibble(rownames = NA)

#select gene information (such as ontology of each gene, chromosome placement, etc)
gene_details=data[["GSE112676_series_matrix.txt.gz"]]@featureData@data

#write clinical and gene data to computer
write_tsv(clindata%>% rownames_to_column(), file = "data/raw/clindata.tsv")
write_tsv(gene_counts2%>% rownames_to_column(), file = "data/raw/gene_counts.tsv")
write_tsv(gene_details%>% rownames_to_column(), file = "data/raw/gene_details.tsv")
