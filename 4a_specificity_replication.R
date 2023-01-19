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
res_pwas <- fread( "Data/Modified/res_pwas.txt")
res_map <- fread("Data/Modified/res_map.txt")
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")
dt_tra <- fread( "Data/Modified/dt_specificity_gene_region.txt")

##################
parameters <- GagnonMR::default_param()
parameters$path <- c(default_param()$path, "/mnt/sda/gagelo01/Projects/Brain_pQTL/Data/Modified/Gtex_vcf/")

##########run analysis and define res_bloodpwmr
options(future.globals.maxSize= 5e9)
plan(multisession, workers = 20, gc = TRUE)
run_all_pqtl_analyses_rosmap_banner_safely <- safely(run_all_pqtl_analyses_rosmap_banner)

arguments <- dt_tra[UniProt %in% c(res_pwas$UniProt, res_map$UniProt),]
arguments[, index := study %>%ifelse(. %in% c("deCODE", "FENLAND"), "blood_protein", .) %>%
            ifelse(. =="eQTLGen", "blood_rna", .) %>%
            ifelse(. == "GTEX", "brain_rna", .)]
arguments[,index_prot:=paste0(index,UniProt)]
arguments <- arguments[, .(UniProt, index_prot, study)]
arguments[,out_wd := ifelse(study == "GTEX", 
                            "/mnt/sda/gagelo01/Projects/Brain_pQTL/Data/Modified/Gtex_vcf/",
                            "/mnt/sda/gagelo01/Vcffile/Server_vcf/")]

res_all <-  future_map(split(arguments, arguments$index_prot), function(x) {
  run_all_pqtl_analyses_rosmap_banner_safely(uniprot = unique(x$UniProt),
                                             out_wd = unique(x$out_wd),
                                             vcffile_out = "/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-1-1/trait-1-1.vcf.gz",
                                             method_list = list("get_uni_cis", "get_coloc"),
                                             study_toinclude = x$study,
                                             dt_gene_region = dt_tra,
                                             parameters = parameters)}, .options = furrr_options(seed = TRUE))

saveRDS(res_all, "Data/Modified/res_specificity.rds")
index <- sapply(res_all, function(x) !is.null(x$result))
res_all <- lapply(res_all[index], function(x) x$result) %>% rbindlist(., fill = TRUE)
res_all <- merge(res_all, dt_tra[,.(id, UniProt, trait, gene_region, hgnc)], by.x = "id.exposure", by.y = "id")

fwrite(res_all, "Data/Modified/res_specificity.txt")

message("This script finished without errors")