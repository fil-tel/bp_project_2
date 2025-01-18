library(mitor)
library(vcfR)
library(tidyverse)
library(grantham)

setwd("UCPH/Year2/block2/analyses/vcf file/")
vcf <- read.vcfR("gnomad.genomes.v3.1.sites.chrM.vcf")

# convert vcf file to dataframe
vcf_df <- data.frame(vcf@fix)



# extract from the info field the allele frequency of the variant (homoplasmic) and the vep field
# those will be used as paramters for filtering
infos <- extract_info_tidy(vcf, info_fields = c("AF_hom", "vep"))

vep_df <- infos$vep %>% str_split(pattern = "\\|")

vep_df <- as.data.frame(do.call(rbind, vep_df))

# header of vep format with | as separator
vep_format <- "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info"
vep_header <- str_split_1(vep_format, "\\|")
# make it double for variants that have an influence on two genes, so from 45 fields we will obtain 90
vep_header <- c(vep_header, paste0(vep_header, "_", 2))

colnames(vep_df) <- vep_header

# delete info column
vcf_df$INFO <- NULL
# bind other information used for filtering
vcf_df <- cbind(vcf_df, AF_hom=infos$AF_hom, vep_df)

vcf_df$POS <- as.numeric(vcf_df$POS)
vcf_df$AF_hom <- as.numeric(vcf_df$AF_hom)

# extract variations corresponding to only coding region for 3D mapping
# I will use the genebank reference that comes with the mitor package

# select only cds
cds_coord <- rCRS_genes_df %>% filter(!startsWith(gene, "T"), !startsWith(gene, "R"))

# PROBABY I DON'T NEED THIS, ENOUGH TO PICK THE MISSENSE VARIANTS
# position to use to "pick" the genes variants
# make_picks <- function(x) cds_coord$start[x]:cds_coord$end[x]
# picks <- unlist(sapply(1:13,make_picks))

# first pick positions, then only missense variants
vcf_df_miss <- vcf_df %>% filter(Consequence == "missense_variant")

# let's have a look at the data
hist(vcf_df_miss$AF_hom, breaks = 2000)

# count the remaining
sum(vcf_df_miss$AF_hom>0.005)

# filter variants with maf>0.001

vcf_df_maf <- vcf_df_miss %>% filter(AF_hom > 0.01)

# let's have a look at the data
hist(vcf_df_maf$AF_hom, breaks = 50)

# convert them to 3 letter code
a_sub <- strsplit(vcf_df_maf$Amino_acids, split = "\\/")
aaa_sub <- lapply(a_sub, seqinr::aaa)

# grantham distance of substitution

gr_score <- lapply(aaa_sub, function(x) grantham_distance(x[1], x[2]))
gr_score_df <- as.data.frame(matrix(unlist(gr_score), ncol = 3, byrow = 3))
colnames(gr_score_df) <- c("REF", "ALT", "Gr_score")

# create df for each prtoein with variable position in protein, amino acid differences, grantham distance
# this will be then used in pymol

df2save <- cbind(vcf_df_maf[c( "SYMBOL", "Protein_position", "AF_hom")], gr_score_df)

write_csv(df2save, "mt_var_01")

