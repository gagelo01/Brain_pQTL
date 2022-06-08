#!/usr/bin/env Rscript
suppressWarnings(suppressPackageStartupMessages({
  library(gwasvcf)
  library(VariantAnnotation)
  library(data.table)
  library(tidyverse)
  library(magrittr)
  library(GagnonMR)
  library(gwasglue)
  library(furrr)
}))
set_bcftools()
setwd("/mnt/sde/gagelo01/Projects/Brain_pQTL")
df_index <- fread("/mnt/sdf/gagelo01/Vcffile/server_gwas_id.txt")

res_pwas <- fread( "Data/Modified/res_pwas.txt")
res_map <- fread("Data/Modified/res_map.txt")
res_all <- rbindlist(list(res_pwas, res_map), fill = TRUE) %>% distinct(.)
res_all <- res_all[, .(exposure, study)] %>% distinct

vec_id <- vector(mode = "logical", length = nrow(res_all))
for(i in 1:nrow(res_all)) {
vec_id[i]<- df_index[trait == res_all[i,]$exposure & consortium == res_all[i,]$study,]
}
ID_exposure <- vec_id %>% unlist(.)

df_exp <- df_index[id %in% ID_exposure, .(trait, id, consortium)]
df_exp[, id := paste0("/mnt/sdf/gagelo01/Vcffile/Server_vcf/", id, "/", id, ".vcf.gz")]
colnames(df_exp) <- c("Gene", "ID_exposure_file", "study")


ID_outcome <- df_index[pmid == 32193382 & !(trait %in% paste0("PC",11:81)),]$id #I take only the ten first principal component
ID_outcome_file <- paste0("/mnt/sdf/gagelo01/Vcffile/Server_vcf/", ID_outcome, "/", ID_outcome, ".vcf.gz")

df_comb <- tidyr::crossing(df_exp, ID_outcome_file)
list_comb <- split(df_comb, seq(nrow(df_comb)))

options(future.globals.maxSize= 1e10)
plan(multisession, workers = 4)

generegion <- GagnonMR::from_genecard_to_generegion( df_exp$Gene %>% unique, window = 1e5)


run_allpqtl_dh<- function(comb, generegion, df_index) {
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  res <- GagnonMR::run_all_pqtl_analyses(comb$ID_exposure_file, comb$ID_outcome_file, chrompos = generegion[comb$Gene], 
                                  method_list = list("get_uni_cis", "get_coloc"))
  res$study <- comb$study
  return(res)
  }
                      
run_allpqtl_dh_safely <- safely(run_allpqtl_dh)
                      
res_coloc <- future_map(list_comb, function(x) {
  run_allpqtl_dh_safely(comb = x, generegion = generegion)},
  .options = furrr_options(seed = TRUE)) 

saveRDS(res_coloc, "Data/Modified/res_coloc_dh")
message("This script finished without errors")


