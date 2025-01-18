# library(devtools)
# install_github("fil-tel/mitor")
library(mitor)
library(magrittr)
library(dplyr)
library(ggplot2)
library(Biostrings)
library(grantham)

setwd("~/UCPH/Year2/block2/analyses/")

# fetch Neanderthals sequences
nean_query <- "(015400[SLEN]:016600[SLEN]) AND txid63221[orgn] AND mitochondrion[FILT]" # 34 results
nean_path <- "data/mtdna/neanderthals"
fetch_seq(nean_query, nean_path, filename = "neanderthals")

# fetch Denisovans sequences
deni_query <- "(015400[SLEN]:016600[SLEN]) AND txid741158[orgn] AND mitochondrion[FILT]" # 8 results
deni_path <- "data/mtdna/denisovans"
fetch_seq(deni_query, deni_path, filename = "denisovans")

# fetch Sima de los Huesos sequences
sima_query <- "(015400[SLEN]:016600[SLEN]) AND txid1425170[orgn] AND mitochondrion[FILT]" # 2 results (repicate)
sima_path <- "data/mtdna/sima"
# n=1 because the two results are the same
fetch_seq(sima_query, sima_path, n=1, filename = "sima")

# fetch all homo sapiens sequences
# I had to modify the query because they changed something in the db and I was fetching neanderthals and denisovans as well
hs_query <- "(015400[SLEN]:016600[SLEN]) AND mitochondrion[FILT] AND txid9606[orgn] NOT txid63221[orgn] NOT txid741158[orgn] NOT NC_012920.1[accn]" 
hs_path <- "data/mtdna/modern_humans"
# fetch all
fetch_seq(hs_query, hs_path, filename = "modern_humans")

# fetch all chimpz
chimp_query <- "(015400[SLEN]:016600[SLEN]) AND txid9598[orgn] AND mitochondrion[FILT])"
chimp_path <- "data/mtdna/chimpz"
# fetch only 10 sequences
fetch_seq(chimp_query, chimp_path, n=10, filename = "chimpz")

# function to clean the record ids in the files and keep only the accession ID (maybe to add in fetch function?)
# Example.
# FROM NC_013993.1 Homo sp. Altai mitochondrion, complete genome TO NC_013993.1

clean_fasta_ids <- function(path2file){
  dnaset <- Biostrings::readDNAStringSet(path2file)
  names(dnaset) <- sapply(strsplit(names(dnaset), split = " "), "[[", 1)
  Biostrings::writeXStringSet(dnaset, filepath = path2file)
}

# clean all
clean_fasta_ids("data/mtdna/denisovans/denisovans.fa")
clean_fasta_ids("data/mtdna/sima/sima.fa")
clean_fasta_ids("data/mtdna/neanderthals/neanderthals.fa")
clean_fasta_ids("data/mtdna/chimpz/chimpz.fa")
clean_fasta_ids("data/mtdna/modern_humans/modern_humans.fa")

# filter out modern humans sequences with IUPAC ambiguities ("[URYSWKMBDHV]") that are not accepted by haplogrep3
# and rewrite the file

mod_hum_set <- Biostrings::readDNAStringSet(paste0(hs_path, "/modern_humans.fa"))

# length(mod_hum_set)
# [1] 67463
mod_hum_set <- mod_hum_set[!grepl("[URYSWKMBDHV]", mod_hum_set)]
# length(mod_hum_set)
# [1] 63397
# ~4000 sequences removed
Biostrings::writeXStringSet(mod_hum_set, "data/mtdna/modern_humans/modern_humans_clean.fa")




# RUN FROM HERE 

mod_hum_set <- Biostrings::readDNAStringSet("data/mtdna/modern_humans/modern_humans_clean.fa") 

# classify modern humans sequences haplogroups using haplogrep3
# maybe better straight in the terminal?

# system("haplogrep3 classify --in data/mtdna/modern_humans/modern_humans_clean.fa --out data/haplogroups_mh --tree phylotree-rcrs@17.2")

# read haplogroups table, in order to sample sequences belonging to each haplogroups from it

haplogroups_df <- read.table("data/haplogroups_mh", header = T)

# filter out low quality sequences (e.g. <85% (idk, just a guess))

# dim(haplogroups_df)
# [1] 63397     5

haplogroups_df <- haplogroups_df %>% filter(Quality>0.85)

# dim(haplogroups_df)
# [1] 58797     5
# ~5000 seqs out

# simplify haplogroups classification considering only higher haplogroups (first letter)

haplogroups_df$Haplogroup <- strtrim(haplogroups_df$Haplogroup, 1)

# let's have a look at the distribution of haplogroups in the database
counts_df <- haplogroups_df %>% group_by(Haplogroup) %>% count()
ggplot(data = counts_df, mapping = aes(x=Haplogroup, y=n)) + geom_histogram(stat = "identity")

# sample 50 sequences (if available !!!) from each haplogroup

hapgroups_sampler <- function(hap_df, n=50){
  ids <- c()
  for (let in LETTERS) {
    extractor <- hap_df$Haplogroup==let
    c <- n
    if(sum(extractor)<n) c <- sum(extractor)
    ids <- c(ids, sample(hap_df$SampleID[extractor], size = c))
  }
  ids
}

sampled_ids <- hapgroups_sampler(haplogroups_df)


# Create Stringset of sequences to align
# start with the revised Cambridge Reference Sequence

seqs2align <- Biostrings::DNAStringSet(rCRS)

# read Denisovans, Neanderthals, Sima and chimpz
deni_set <- Biostrings::readDNAStringSet("data/mtdna/denisovans/denisovans.fa")
sima_set <- Biostrings::readDNAStringSet("data/mtdna/sima/sima.fa")
nean_set <- Biostrings::readDNAStringSet("data/mtdna/neanderthals/neanderthals.fa")
chimpz_set <- Biostrings::readDNAStringSet("data/mtdna/chimpz/chimpz.fa")

# add all sequences, together with sampled modern humans

seqs2align <- c(seqs2align, mod_hum_set[sampled_ids], nean_set, deni_set,  sima_set, chimpz_set)

# create directory and save the sequence there

dir.create("data/mtdna/seq2align", recursive = TRUE)
Biostrings::writeXStringSet(seqs2align, "data/mtdna/seq2align/seq2align.fa")


# NOW RUN MSA, I'LL DO NOT IN R CAUSE FASTER AND MORE STABLE

# clustalo --in seq2align.fa --out msa.fa --output-order input-order --threads 2

my_msa <- Biostrings::readDNAStringSet("data/mtdna/seq2align/msa.fa")

# get names of sequences

deni_ids <- names(Biostrings::readDNAStringSet("data/mtdna/denisovans/denisovans.fa"))
sima_ids <- names(Biostrings::readDNAStringSet("data/mtdna/sima/sima.fa"))
nean_ids <- names(Biostrings::readDNAStringSet("data/mtdna/neanderthals/neanderthals.fa"))
chimpz_ids <- names(Biostrings::readDNAStringSet("data/mtdna/chimpz/chimpz.fa"))
mh_ids <- setdiff(names(my_msa), c(deni_ids, sima_ids, nean_ids, chimpz_ids))

# modern humans to do later

adj_coord_df <- genes_coord(my_msa)



# helper function to get the ids of the sequences to remove

get_stop_ids <- function(prot_msa){
  ids <- c()
  # matrix of the variale positions
  var_mat <- find_var_pos(prot_msa)
  ids <- apply(var_mat, 2, function(x){
    if("*" %in% x){
      c(ids, names(x[x=="*"]))
    }
  })
  unique(ids)
}

get_ids <- function(cds_msa){
  ids <- c()
  # matrix of the variale positions
  var_mat <- find_var_pos(cds_msa, type = "DNA")
  ids <- apply(var_mat, 2, function(x){
    if(x["NC_012920"]=="-"){
      c(ids, names(x[x!="-"]))
    } else if("-" %in% unique(x)){
      c(ids, names(x[x=="-"]))
    }
  })
  unique(ids)
}

# function to filter out seqeunces with frameshift in cds
# this because I am only interested in non-synonimous mutations
# it returns the MSA without the "incriminated" sequences

clean_msa <- function(msa){
  adj_coord_df <- genes_coord(msa)
  cds_list <- extract_cds(my_msa, adj_coord_df)
  ids2rm <- c()
  for (name in names(cds_list)) {
    ids2rm <- c(ids2rm, get_ids(cds_list[[name]]))
  }
  ids2rm <- unique(unlist(ids2rm))
  # print(ids2rm)
  msa_filt1 <- msa[setdiff(names(msa), ids2rm)]
  # clean for premature stop codons
  
  cds_list <- extract_cds(msa_filt1, adj_coord_df)
  protein_list <- lapply(cds_list, translation)
  
  ids2rm <- c()
  for (name in names(protein_list)) {
    ids2rm <- c(ids2rm, get_stop_ids(protein_list[[name]]))
  }
  ids2rm <- unique(unlist(ids2rm))
  # print(ids2rm)
  msa_filt1[setdiff(names(msa_filt1), ids2rm)]
}

# length(my_msa)
# 1288

my_msa_filt <- clean_msa(my_msa) 

# length(my_msa_filt)
# 1271 >>> 17 removed

# update names

tot_ids <- names(my_msa_filt)
deni_ids <- intersect(tot_ids, deni_ids)
sima_ids <- intersect(tot_ids, sima_ids)
nean_ids <- intersect(tot_ids, nean_ids)
chimpz_ids <- intersect(tot_ids, chimpz_ids)
mh_ids <- intersect(tot_ids, mh_ids)


cds_list <- extract_cds(my_msa_filt, adj_coord_df)

protein_list <- lapply(cds_list, translation)


# Function to include in package 

find_variants_AA <- function(msa, target, ref="NC_012920"){
  # ref seq goes in front
  if(!target %in% names(msa)) stop("The target sequence is not present in the MSA. Check the argument you passed.")
  if(!ref %in% names(msa)) stop("The reference sequence is not present in the MSA. Check the argument you passed.")
  seqs <- as.matrix(msa[c(ref,target)])
  variants_list <- apply(seqs, 2, unique)
  gaps <- unlist(gregexpr("-",as.character(msa[[ref]])))
  if(gaps==-1) gaps <- NULL
  mutations <- c()
  i <- 1
  # variants_list
  while (i<=length(variants_list)) {
    # if in the position there is only 1 character (no difference) skip the iteration
    if(length(variants_list[[i]])==1 |  "X" %in% variants_list[[i]]){
      i <- i+1
      next
    }
    # Insertion case: if there is a gap in the reference
    else if(variants_list[[i]][1]=="-"){
      # position after the beginning of the gap
      if(i+1<=length(variants_list)){
        c <- i+1
      }
      # if the insertion is longer than one base we need to check how long it is
      while (variants_list[[c]][1]=="-" & c+1<=length(variants_list) & length(variants_list[[c]])==2) {
        c <- c+1
      }
      # adjust the insertion indices according to the position in the ref seq not in the msa, as always
      ins_start <- i-sum(gaps<i)-1
      ins_end <- ins_start+1
      # extract the bases inserted in the target seq
      ins_bases <- paste0(sapply(variants_list[i:(c-1)], function(x) x[2]), collapse = "")
      mutations <- c(mutations, paste0(ins_start, "_", ins_end, "ins", ins_bases))
      # this is needed in case we have a gap in the last position, otherwise if don't do this we get stuck in the loop
      if(!is.na(variants_list[[c]][2]) & variants_list[[c]][2]=="-" & c+1>length(variants_list)){
        i <- c+1
      }
      else{
        i <- c
      }
    }
    # Deletion case: if there is a gap in the target
    else if(variants_list[[i]][2]=="-"){
      # position after the beginning of the gap
      if(i+1<=length(variants_list)){
        c <- i+1
      }
      # if the deletion is longer than one base we need to check how long it is
      # the other conditions are to avoid to go out of bound when if at the end of the alignment
      while (!is.na(variants_list[[c]][2]) & variants_list[[c]][2]=="-" & c+1<=length(variants_list) & length(variants_list[[c]])==2) {
        c <- c+1
      }
      # adjust the insertion indices according to the position in the ref seq not in the msa, as always
      del_start <- i-sum(gaps<i)
      del_end <- del_start+(c-i)-1
      # one base del
      if(del_start==del_end){
        del_bases <- variants_list[[i]][1]
        mutations <- c(mutations, paste0(del_start, "del", del_bases))
      }
      # more than one base deletion
      else{
        # extract the bases inserted in the target seq
        del_bases <- paste0(sapply(variants_list[i:(c-1)], function(x) x[1]), collapse = "")
        mutations <- c(mutations, paste0(del_start, "_", del_end, "del", del_bases))
      }
      # this is needed in case we have a gap in the last position, otherwise if don't do this we get stuck in the loop
      if(!is.na(variants_list[[c]][2]) & variants_list[[c]][2]=="-" & c+1>length(variants_list)){
        i <- c+1
      }
      else{
        i <- c
      }
    }
    # subsituition case
    else{
      # adjust index according to the msa
      sub_i <- i-sum(gaps<i)
      mutations <- c(mutations, paste0(sub_i, variants_list[[i]][1], ">", variants_list[[i]][2]))
      i <- i+1
    }
  }
  mutations
}


# find consensus sequences for each group

cons_mh <- lapply(protein_list, function(x) Biostrings::consensusString(x[mh_ids]))
cons_nean <- lapply(protein_list, function(x) Biostrings::consensusString(x[nean_ids]))
cons_deni <- lapply(protein_list, function(x) Biostrings::consensusString(x[deni_ids]))
cons_chimpz <- lapply(protein_list, function(x) Biostrings::consensusString(x[chimpz_ids]))
cons_sima <- lapply(protein_list, function(x) Biostrings::consensusString(x[sima_ids]))

# now I want to detect for each linage, which amino acid are derived in that specific
# branch of the tree, that differs from the ancestral state (chimp)

detect_derived_aa <- function(msa, prot_name = NULL){
  if(is.null(prot_name)) stop("Provide a protein name!")
  var_pos <- find_var_pos(msa)
  
  
  
}


# test

msa_ex <- Biostrings::AAStringSet(c("Chimp"=cons_chimpz$ND1, "Humans"=cons_mh$ND1, "Neanderhals"=cons_nean$ND1, "Denisovans"=cons_deni$ND1, "Sima"=cons_sima$ND1))

var_pos <- find_var_pos(msa_ex)

spot <- function(col){
  if(length(unique(col[names(col)!="Chimp"]))==1 & !("X" %in% col[names(col)!="Chimp"]))
    return("Homo")
  else if(length(unique(col[names(col)!="Chimp"]))==1)
}

apply(var_pos, 2, function(x) spot(x))
