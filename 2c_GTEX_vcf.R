#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(GagnonMR)
library(furrr)
setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL")

#load objects
dt_res <- fread( "Data/Modified/dt_res_pwas.txt")
res_pwas <- fread( "Data/Modified/res_pwas.txt")
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
res_map <- fread( "Data/Modified/res_map.txt")
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")

#determine the list of gene and the tissue

uniprot_vec <- unique(c(dt_res[!is.na(pval.exposure.wald), unique(UniProt)], res_map$UniProt))
vec_gene <- dt_gene_region[UniProt %in% uniprot_vec, ]$hgnc %>% unique
vec_tissue <-"Brain_Frontal_Cortex_BA9"
plan(multicore, workers = 9, gc = TRUE) 
plan(sequential)

#NO need to change from here
options(future.globals.maxSize= 5e9)
get_eQTL_safely <- safely(get_eQTL)


df_eqtl <-  future_map(as.list(vec_gene), function(x) {
  GagnonMR::get_eQTL(tissue = vec_tissue, gene = x, mywindow = 2e5)  }, .options = furrr_options(seed = TRUE)) %>% rbindlist(.,fill = TRUE)


fwrite(df_eqtl, "Data/Modified/df_eqtl.txt")

########format gtex and make it build 37
translation <- fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b38_rsid_maf_small.txt")
traduction <-  fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt")
traduction[, EUR := EUR %>% ifelse(.==0,0.001,. ) %>% ifelse(.==1, 0.999, .)]
traduction[, maf := NULL]
gencode <- fread("/home/couchr02/Mendel_Commun/Christian/GTEx_v8/gencode.v26.GRCh38.genes.txt")

exposures_formated <- format_gtex_data(exposures = df_eqtl, translation = translation, traduction = traduction, gencode = gencode)
exposures_formated<- separate(exposures_formated, col = "id.exposure", into = c("id.exposure", "gene.exposure"), sep = "-")
exposures_gtex <- exposures_formated
fwrite(exposures_gtex, "Data/Modified/exposures_gtex_hyprcoloc.txt")


##############make gtex as vcf###################
info <- data.table::data.table(trait =exposures_gtex[, unique(exposure)])

newrow <- data.table(id = paste0("eqtl-74-", 1:info[,.N]), trait = info$trait, group_name = "public",
                     year = 2020, author = "The GTEX Consortium",consortium = "GTEX",
                     sex = "Males and Females", population = "European",
                     initial_build = "HG38/GRCh38", unit ="SD", nsnp = "dummy", sample_size = 10708,
                     category = "protein", pmid = 34648354, ncase = NA,sd = 1, note = "dummy", ncontrol = NA,
                     adjustments = "dummy")

df_index <- df_index[!grepl("eqtl-74-", id), ]
df_index <- rbind(df_index, newrow)

exposures_gtex <- exposures_gtex[!(is.na(se.exposure)|is.na(eaf.exposure)),]
#eqtl-4
newrow_split <- split(newrow, newrow$id)
dir.create("Data/Modified/Gtex_vcf")
map(newrow_split, function(x) { 
  GagnonMR::formattovcf_createindex(all_out = exposures_gtex[exposure==x$trait,],
                                    snp_col = "SNP",
                                    outcome_name = x$trait,
                                    beta_col = "beta.exposure",
                                    se_col = "se.exposure",
                                    pval_col = "pval.exposure",
                                    eaf_col = "eaf.exposure",
                                    effect_allele_col = "effect_allele.exposure",
                                    other_allele_col =  "other_allele.exposure",
                                    ncase_col = NULL,
                                    ncontrol_col = NULL,
                                    samplesize_col = "samplesize.exposure",
                                    chr_col = "chr.exposure",
                                    pos_col = "pos.exposure",
                                    units = "SD",
                                    traduction = traduction, #already harmonised with traduction
                                    out_wd = "/mnt/sda/gagelo01/Projects/Brain_pQTL/Data/Modified/Gtex_vcf",
                                    df_index = df_index,
                                    group_name = "public",
                                    year = 2020,
                                    author = "The GTEX Consortium",
                                    consortium = "GTEX",
                                    sex = "Males and Females",
                                    population = "European",
                                    initial_build = "HG38/GRCh38",
                                    category = "eqtl",
                                    pmid = 32913098,
                                    note = "",
                                    should_create_id = FALSE,
                                    ID = x$id) 
})
message("This script finished without errors")


fwrite(newrow, "Data/Modified/newrow_eqtlgtex.txt")
