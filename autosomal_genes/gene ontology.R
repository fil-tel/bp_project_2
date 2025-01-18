library(readxl)
library(gprofiler2)
library(biomaRt)
library(dplyr)
setwd("/home/filippo/UCPH/Year2/block2/analyses/autosomal_genes/")

amhir <- read_xlsx("AMHIRgenes.xlsx", col_names = FALSE)
colnames(amhir) <- "Genes"
# write(amhir$Genes, "amhir_genes")

amhir_vec <- amhir$Genes

# GO enrichment analysis
gostres <- gost(query = amhir_vec, organism = "hsapiens", sources = "GO")
head(gostres$result)

p <- gostplot(gostres, capped = FALSE, interactive = TRUE)
p

# Retrieve all genes associated with a go term and his child terms
# I want to retrieve all the genes associated with mitochondria: GO:0005739

ensembl <- useEnsembl(biomart="ensembl", dataset='hsapiens_gene_ensembl')
# get all the genes as a dataframe
mtgo_genes <- getBM(attributes = c('ensembl_gene_id','external_gene_name', 'name_1006'), 
      filters = 'go', 
      values = 'GO:0005739', 
      mart = ensembl)

# look if genes in my list match with those genes
my_mt_genes <- c()

# select unique because gene names are repeated
unique_mtgo_genes <- unique(mtgo_genes$external_gene_name)
for(gene in amhir_vec){
  if(gene %in% unique_mtgo_genes){
    my_mt_genes <- c(my_mt_genes, gene)
  }
    
}

my_mt_genes

# filter the dataframe according to those genes

amhir_mt_genes <- mtgo_genes %>% filter(external_gene_name %in% my_mt_genes)

# now I try to run a GOEA on those mt genes using as background all the mt genes
# to see if there is an enrichment of certain mt function in the mt introgressed genes

gostres_mt <- gost(query = my_mt_genes, organism = "hsapiens", sources = "GO", custom_bg = unique_mtgo_genes)
head(gostres$result)

p <- gostplot(gostres, capped = FALSE, interactive = TRUE)
p

# 
# 
# 
# 
# write list
write(unique(amhir_mt_genes$ensembl_gene_id), "ensembl_list_mt")
# write mt background
write(unique(mtgo_genes$ensembl_gene_id), "mt_background")
