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
outcome_file <- fread("/mnt/sdf/gagelo01/Dietary_habits/outcome_file-translation.csv")
colnames(outcome_file) <- c("file", "outcome")
defiles <- list.files("/mnt/sdf/gagelo01/Dietary_habits/")
defiles <- defiles[grepl(".gz",defiles)]
mmm<- sub("BOLTlmm_UKB_genoQCEURN455146_v3_diet_", "", defiles) %>% sub("_BgenSnps_mac20_maf0.005_info0.6.gz", "", .)
outcome_file <- outcome_file[file %in% mmm, ] 
outcome_file[, outcome := gsub(" |-", "_", outcome) %>% gsub("\\.|,|:|", "", .) %>% gsub("\\+", "and", .) %>% gsub("/", "_", .) %>% gsub("\\(", "", .) %>% gsub("\\)", "", .)]

newrow <- data.frame(id = NA, trait = outcome_file$outcome, group_name = "Dietary habits", year = 2021, author = "B. Cole, Joanne",
                     consortium = "UKBbiobank",
                     sex = "Males and Females", population = "European", unit = "SD", nsnp = NA, sample_size = 449210,
                     initial_build = "HG19/GRCh37", category = "Trait", pmid = 32193382, ncase = NA,
                     sd = 1, note = NA, ncontrol = NA)



df_index <- rbindlist(list(df_index, newrow), fill = TRUE)
df_index[is.na(id), id := paste0("trait-3-", 1:nrow(newrow))]
dt_index <- df_index

traduction <-  fread("/mnt/sde/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt")
traduction[ , EUR := EUR %>% ifelse(. == 0, 0.001, .) %>% ifelse(. == 1, 0.999, .)]

#
wrapper_formatvcf <- function(ID, outcome_file, traduction, dt_index) {
file <- paste0("BOLTlmm_UKB_genoQCEURN455146_v3_diet_", outcome_file[outcome == df_index[id == ID, ]$trait,]$file, "_BgenSnps_mac20_maf0.005_info0.6.gz")
all_out <- fread(paste0("/mnt/sdf/gagelo01/Dietary_habits/", file))
all_out <- all_out[grepl("^rs", SNP),]

formattovcf_createindex(all_out = all_out,
                        snp_col = "SNP",
                        outcome_name = outcome_file[outcome == df_index[id == ID, ]$trait,]$outcome,
                        beta_col = "BETA",
                        se_col = "SE",
                        pval_col = "P_BOLT_LMM",
                        eaf_col = "A1FREQ",
                        effect_allele_col = "ALLELE1",
                        other_allele_col = "ALLELE0",
                        ncase_col = NULL,
                        ncontrol_col = NULL,
                        samplesize_col = 449210,
                        chr_col = "CHR",
                        pos_col = "BP",
                        units = "SD",  
                        traduction = traduction,
                        out_wd = "/mnt/sdf/gagelo01/Vcffile/Server_vcf",
                        df_index = dt_index,
                        group_name = "Dietary habits",
                        year = 2021, 
                        author = "B. Cole, Joanne",
                        consortium = "UKBbiobank",
                        sex = "Males and Females",
                        population = "European", 
                        initial_build = "HG19/GRCh37", 
                        category = "Trait", 
                        pmid = 32193382,
                        note = NA,
                        should_create_id = FALSE,
                        ID = ID)
}



options(future.globals.maxSize= 5e9)
plan(multisession, workers = 5)

vec_id <-df_index[trait %in% outcome_file$outcome,]$id

outcome_file_copy <- outcome_file
future_map(as.list(vec_id), function(x) {
  wrapper_formatvcf(ID = x, outcome_file = outcome_file_copy, traduction = traduction, dt_index = df_index)},
 .options = furrr_options(seed = TRUE))


fwrite(dt_index, "/mnt/sdf/gagelo01/Vcffile/server_gwas_id.txt")
message("This script finished without errors")
  