#!/usr/bin/env Rscript
#run_analysis
library(data.table)
library(tidyverse)
library(GagnonMR)
setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL")

df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")
newrow <- fread( "Data/Modified/newrow_eqtlgtex.txt")
newrow <- newrow[,.(id, trait)]
newrow <- separate(newrow, "trait", into = c("tissue", "hgnc_symbol"), sep = "-", remove = FALSE)

#####Fenland Decode
dt_fd <- df_index[grepl("^prot-7|^prot-5", id), .(id, trait, consortium, note)]
dt_fd[, trait := gsub("^fenland", "", trait)]
setnames(dt_fd, c("consortium", "note"), c("study", "UniProt"))
dt_fd <- merge(dt_fd, dt_gene_region[, .SD[1], by = "UniProt"][,.(gene_region, hgnc, UniProt)], by = "UniProt")

##########eQTLGen#########
mart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
trans <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id",
                        values = df_index[grepl("^eqtl-3-", id), trait], mart = mart)

setDT(trans)
tra <- trans
tra <- merge(tra, distinct(dt_gene_region[, .(hgnc, UniProt, gene_region)]), by.x = "hgnc_symbol", by.y = "hgnc")
tra <- merge(tra, df_index[grepl("^eqtl-3-", id), .(id, trait, consortium)], by.x = "ensembl_gene_id", by.y = "trait")
setnames(tra, c("consortium", "ensembl_gene_id", "hgnc_symbol"), c("study", "trait", "hgnc"))

#####GTEX#######
gt<-newrow
setDT(gt)
gt <- merge(gt, distinct(dt_gene_region[, .(hgnc, UniProt, gene_region)]), by.x = "hgnc_symbol", by.y = "hgnc")
gt[,consortium := "GTEX"]
gt[,tissue:= NULL]
setnames(gt, c("consortium",  "hgnc_symbol"), c("study", "hgnc"))

########create dt_tra  #####
dt_tra <- rbindlist(list(dt_fd, tra, gt), fill = TRUE)
fwrite(dt_tra, "Data/Modified/dt_specificity_gene_region.txt")

message("This script finished without errors")
