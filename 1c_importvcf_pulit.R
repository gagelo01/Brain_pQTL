#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)

gwasvcf::set_bcftools()
gwasvcf::set_plink()

setwd("/mnt/sde/gagelo01/Projects/Brain_pQTL")

sumstat <- fread("/mnt/sde/gagelo01/Projects/Brain_pQTL/Data/Raw/Outcome/bmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz")
sumstat[,SNP := gsub(":.", "", SNP)]

traduction <- fread("/mnt/sde/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt")
traduction[, EUR := EUR %>% ifelse(.==0,0.001,. ) %>% ifelse(.==1, 0.999, .)]
traduction[, maf := NULL]


df_index <- fread("/mnt/sdf/gagelo01/Vcffile/server_gwas_id.txt")


GagnonMR::formattovcf_createindex(all_out = sumstat,
                                  snp_col = "SNP",
                                  outcome_name = "bmi_ukbgiant",
                                  beta_col = "BETA",
                                  se_col = "SE",
                                  pval_col = "P",
                                  eaf_col = "Freq_Tested_Allele",
                                  effect_allele_col = "Tested_Allele",
                                  other_allele_col =  "Other_Allele",
                                  ncase_col = NULL,
                                  ncontrol_col = NULL,
                                  samplesize_col = "N",
                                  chr_col = "CHR",
                                  pos_col = "POS",
                                  units = "SD",
                                  traduction = traduction,
                                  out_wd = "/mnt/sdf/gagelo01/Vcffile/Server_vcf",
                                  df_index = df_index,
                                  group_name = "",
                                  year = 2018,
                                  author = "Pulit Sara L",
                                  consortium = "GIANT and UKBiobank",
                                  sex = "Males and Females",
                                  population = "European",
                                  initial_build = "HG19/GRCh37",
                                  category = "Trait",
                                  pmid = 30239722,
                                  note = "",
                                  should_create_id = FALSE,
                                  ID = "trait-1-1")
