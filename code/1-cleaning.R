# load libraries
library(tidyverse)
  


#load raw gene expression data
raw_gene_counts=read_tsv(file="data/raw/gene_counts.tsv")

#add rownames to gene counts
raw_gene_counts=raw_gene_counts%>% 
  remove_rownames %>% 
  column_to_rownames(var="rowname")
  
#convert gene counts to matrix
raw_gene_counts<- as.matrix(raw_gene_counts)

#transpose and convert to tibble
raw_gene_counts=t(raw_gene_counts)%>%
  as_tibble(rownames=NA)%>%
  rownames_to_column("rowname")

#load clinical data 
clindata=read_tsv(file="data/raw/clindata.tsv")

#clean clinical data
#rename to simpler names and convert diagnosis to binary
clindata=clindata%>%
  rename_at(vars(contains(":ch1")), ~ str_remove(., ":ch1*"))%>%
  select(c("diagnosis",'rowname'))%>%
  mutate(diagnosis = recode(diagnosis, 
                    "CON" = 0, 
                    "ALS" = 1))

#join clinical data to gene data
ALS_data=merge(clindata,raw_gene_counts, by = "rowname")%>%
  column_to_rownames(var="rowname")

#select only genes with no NA values for functional and pathway associated information
#also use regex to clean functional and pathway associated information text
gene_info_clean=read_tsv("data/raw/gene_details.tsv")%>%
  drop_na(c("Ontology_Component","Ontology_Process","Ontology_Function"))%>%
  mutate(Ontology_Component = str_remove(Ontology_Component, '\\[.*'))%>%
  mutate(Ontology_Process = str_remove(Ontology_Process, '\\[.*'))%>%
  mutate(Ontology_Function = str_remove(Ontology_Function, '\\[.*'))%>%
  mutate(Definition = str_remove(Definition, 'Homo sapiens*'))%>%
  mutate(Definition = str_remove(Definition, ',.*'))%>%
  rename(c("Feature ID"="ID","Associated Cellular Component"="Ontology_Component","Associated Cellular Function"="Ontology_Function",
           "Associated Cellular Process"="Ontology_Process","Extended Gene Name"="Definition","Gene Acronym"="ILMN_Gene"))%>%
  select(c('Feature ID',"Gene Acronym","Extended Gene Name",'Associated Cellular Component','Associated Cellular Process',"Associated Cellular Function"))

clean_genes=gene_info_clean%>%
  pull(`Feature ID`)


#select gene expression data of genes with associated data
ALS_data=ALS_data%>%
  select(c('diagnosis',all_of(clean_genes)))

#calculate variance of each gene
variances <- apply(X=ALS_data%>%select(-diagnosis), MARGIN=2, FUN=var)

# sort genes by variance, get top 2500
sorted_variance <- sort(variances, decreasing=TRUE, index.return=TRUE)$ix[1:2501]

ALS_data_final=ALS_data[c(1,sorted_variance)]%>%
  select(-diagnosis.1)

#get top 2500 gene info
gene_info_clean=gene_info_clean%>%
  filter(`Feature ID` %in% colnames(ALS_data_final))

#save cleaned data to file
write_tsv(ALS_data_final, file = "data/clean/ALS_data.tsv")
write_tsv(gene_info_clean, file = "data/clean/gene_info.tsv")

