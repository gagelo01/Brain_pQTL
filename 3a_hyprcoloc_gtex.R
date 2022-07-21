#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(GagnonMR)
library(furrr)
setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL")

#load objects
# gencode <- fread("/home/couchr02/Mendel_Commun/Christian/GTEx_v8/gencode.v26.GRCh38.genes.txt")
# tissue_dir <- list.files("/home/couchr02/Mendel_Commun/Nicolas/GTEx_V8/GTEx_EUR_Analysis_v8_eQTL_expression_matrices")
# tissue_loop <- sub(".v8.normalized_expression_EUR_chr1_22.bed","",tissue_dir)
res_pwas <- fread( "Data/Modified/res_pwas.txt")
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
res_map <- fread( "Data/Modified/res_map.txt")
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")

#determine the list of gene and the tissue
uniprot_vec <- unique(c(res_map$UniProt, res_pwas$UniProt)) 
vec_gene <- dt_gene_region[UniProt %in% uniprot_vec, ]$hgnc %>% unique
vec_tissue <-"Brain_Frontal_Cortex_BA9"
plan(multicore, workers = 9, gc = TRUE) #
plan(sequential)

#NO need to change from here
options(future.globals.maxSize= 5e9)
get_eQTL_safely <- safely(get_eQTL)


df_eqtl <-  future_map(as.list(vec_gene), function(x) {
  GagnonMR::get_eQTL(tissue = vec_tissue, gene = x, mywindow = 2e5)  }, .options = furrr_options(seed = TRUE)) %>% rbindlist(.,fill = TRUE)


fwrite(df_eqtl, "Data/Modified/df_eqtl.txt")

########format gtex and make it build 38
translation <- fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b38_rsid_maf_small.txt")
traduction <-  fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt")
gencode <- fread("/home/couchr02/Mendel_Commun/Christian/GTEx_v8/gencode.v26.GRCh38.genes.txt")

exposures_formated <- format_gtex_data(exposures = df_eqtl, translation = translation, gencode = gencode)

exposures_formated<-merge(exposures_formated, traduction[,.(rsid, chr, position)], by.x = "SNP", by.y = "rsid")
exposures_formated[, chr.exposure := chr]
exposures_formated[, pos.exposure := position]
exposures_formated[, chr := NULL]
exposures_formated[,position := NULL]
exposures_formated[, gene.exposure := NULL]
exposures_formated<- separate(exposures_formated, col = "id.exposure", into = c("id.exposure", "gene.exposure"), sep = "-")
exposures_gtex <- exposures_formated
fwrite(exposures_gtex, "Data/Modified/exposures_gtex_hyprcoloc.txt")

message("This script finished without errors")
