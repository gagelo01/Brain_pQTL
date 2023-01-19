#!/usr/bin/env Rscript
#run_analysis
library(data.table)
library(tidyverse)
library(GagnonMR)
library(tictoc)
library(furrr)
setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL")
source("Analysis/run_all_pqtl_analyses_rosmap_banner.R")

df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
ldref<- "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
snp_bim<-fread(paste0(ldref, ".bim"), select = 2)$V2

###
somavec <- df_index[pmid == 34239129, trait]
dt_aptamers <- as.data.table(readat::aptamers)
all( somavec %in% dt_aptamers$AptamerId)
dtpath <- dt_aptamers[AptamerId %in% somavec, .(AptamerId, UniProt, EntrezGeneSymbol),]
vecpos <- GagnonMR::from_genecard_to_generegion(dtpath$EntrezGeneSymbol, window = 2e5)
vecpos[sapply(vecpos, function(x) length(x)>1)]<-NA #those that have two different chromosomes transform in NA
vecpos[sapply(vecpos, function(x) grepl(paste(c(letters,LETTERS),collapse = "|"), x))]<-NA
vecpos <- unlist(vecpos)
dtpos <- data.table(gene_name = names(vecpos), region = vecpos)
dtpos <- merge(dtpath, distinct(dtpos), by.x = "EntrezGeneSymbol", by.y = "gene_name")

df_index_copy <- df_index[pmid %in% c(33571421, 34239129), ]
df_index_copy <- separate(df_index_copy,  col = "note", into = c("hgnc", "chrompos"), sep = "_")
dt_gene_region <- df_index_copy[, .(id, trait, chrompos, hgnc, consortium, pmid)]
setnames(dt_gene_region, c("chrompos", "consortium"), c("gene_region", "study"))
dt_gene_region <- merge(dt_gene_region, dtpos, by.x = "trait", by.y = "AptamerId", all.x = TRUE)
dt_gene_region[is.na(gene_region), gene_region := region]
dt_gene_region[pmid == 34239129, study := "yang"]
dt_gene_region[hgnc == "", hgnc := EntrezGeneSymbol]
dt_gene_region[is.na(UniProt), UniProt := trait]
dt_gene_region[,c("pmid", "region", "EntrezGeneSymbol"):=NULL]
####################
dt_gene_region[, gene_region := ifelse(is.na(gene_region), unique(gene_region[!is.na(gene_region)]), gene_region), by = "UniProt"]

 
dt_gene_region <- separate(dt_gene_region, col = "gene_region", into = c("chr", "start", "end"), remove = FALSE)
if(!dt_gene_region[, length(unique(chr)), by = "UniProt"][,any(V1>1)]) {
  dt_gene_region[,newgene_region := paste0(unique(chr),":",min(start), "-",max(end)), by = "UniProt"]
  dt_gene_region[,gene_region:=newgene_region]
  dt_gene_region[,c("chr","start", "end", "newgene_region"):=NULL]
}

#####include ensemble gene id
k <- dt_gene_region[,.(hgnc, id)]
nom_var<-"hgnc"
num<-k[,str_count(get(nom_var), pattern = " ")%>%max+1]
k <- separate(data = k, col = nom_var, sep = " ", into = paste0(nom_var, 1:num))
k <- melt(k, id.vars = c("id"), measure.vars = paste0(nom_var, 1:num))
k<- k[!is.na(value),]
mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
dt_ensg <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "hgnc_symbol", values = unique(k$value), mart = mart)     
setDT(dt_ensg)
dt_ensg <- dt_ensg[, .SD[1] ,by = "hgnc_symbol"] #arbitrarily pick a ENSG synonym for the summary if multiple.
setnames(dt_ensg, c("hgnc_symbol", "ensembl_gene_id"), c("hgnc","ensg"))
k <- merge(dt_ensg, k, by.x = "hgnc", by.y = "value", all.y = TRUE)
k <- dcast(k, id ~ variable, value.var = c("hgnc", "ensg"))
k <- unite(k, col = "hgnc", paste0("hgnc_",nom_var, 1:num), sep = " ")
k[,hgnc := gsub(" NA", "", hgnc)]
k <- unite(k, col = "ensg", paste0("ensg_",nom_var, 1:num), sep = " ")
k[,ensg := gsub(" NA", "", ensg)]
dt_gene_region <- merge(dt_gene_region, k, by = c("id", "hgnc"), all.x = TRUE)

######
dt_gene_region_copy <-dt_gene_region
uniprot <- dt_gene_region$UniProt %>% unique
set.seed(seed = 2312)
uniprot <- sample(x = uniprot, size = length(uniprot), replace = FALSE)
parameters <- GagnonMR::default_param()
parameters$snp_bim<-snp_bim

run_all_pqtl_analyses_rosmap_banner_safely <- safely(run_all_pqtl_analyses_rosmap_banner)

options(future.globals.maxSize= 5e9)
plan(multisession, workers = 30, gc = TRUE) #plan(sequential)

tic()

res_all <-  future_map(as.list(uniprot), function(x) {
  run_all_pqtl_analyses_rosmap_banner_safely(uniprot = x, 
                                             vcffile_out = "/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-1-1/trait-1-1.vcf.gz",
                                             method_list = list("get_uni_cis", "get_coloc"),
                                             dt_gene_region = dt_gene_region_copy,
                                             parameters = parameters)}, .options = furrr_options(seed = TRUE))

saveRDS(res_all, "Data/Modified/safely_res_all5.rds")

toc() 

fwrite(dt_gene_region, "Data/Modified/dt_gene_region.txt")
message("This script finished without errors")

  


