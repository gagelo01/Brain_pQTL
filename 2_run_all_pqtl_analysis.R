#!/usr/bin/env Rscript
#run_analysis
library(data.table)
library(tidyverse)
library(GagnonMR)
library(tictoc)
library(furrr)

setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL")

# k_rosmap <- fread( "Data/Modified/k_rosmap.txt")
# k_banner <- fread("Data/Modified/k_banner.txt")


# ensembl_human <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
#                                   host = "grch37.ensembl.org", path = "/biomart/martservice", 
#                                   dataset = "hsapiens_gene_ensembl")
# 
# gencode_copy <- data.table::fread("/home/couchr02/Mendel_Commun/Nicolas/GTEx/gencode.v19.genes.v7.patched_contigs.txt")
# gencode_copy <- gencode_copy[gene_name != "NAA38", ]
# gencode_copy <- gencode_copy[!(gene_name == "LSP1" & chr == 13), ]
# gencode_copy <- rbind(gencode_copy, list(17, 7759999, 7788738, "ENSG00000128534.3", "NAA38", "protein_coding"))
# gencode_copy <- gencode_copy[!(gene_name %in% gencode_copy[, length(unique(chr)),by = "gene_name"][V1 > 1]$gene_name), ] #remove every with more than one chromosome
# 
# gene_region <- GagnonMR::from_genecard_to_generegion(genecard_name = df_index_copy[hgnc != "nonavailable", unique(hgnc)],
#                                                      window = 2e5, 
#                                                      ensembl = ensembl_human, 
#                                                      gencode = gencode_copy)
# dt_gene_region <- data.table(hgnc = names(gene_region), gene_region = unname(gene_region)) 
# dt_gene_region <- merge(dt_gene_region, df_index_copy[,.(hgnc, trait)], by = "hgnc")
# dt_gene_region <- distinct(dt_gene_region)
# 
# dtna<-df_index_copy[hgnc == "nonavailable" | trait %in% dt_gene_region[is.na(gene_region)], .(trait, chrompos, hgnc, consortium)]
# setnames(dna, c("chrompos", "consortium"), c("gene_region", "study"))
# dt_gene_region <- rbind(dt_gene_region, dtna)
# 
# df_index_copy[trait %in% dt_gene_region[is.na(gene_region), trait], ]
# 
# dt_gene_region[,.N] -   df_index_copy$trait %>% unique %>% length


df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
df_index_copy <- df_index[pmid == 33571421, ]
df_index_copy <- separate(df_index_copy,  col = "note", into = c("hgnc", "chrompos"), sep = "_")
dt_gene_region <- df_index_copy[, .(trait, chrompos, hgnc, consortium)]
setnames(dt_gene_region, c("chrompos", "consortium"), c("gene_region", "study"))
dt_gene_region_copy <-dt_gene_region
run_all_pqtl_analyses_rosmap_banner <- function(ID, 
                                      out_wd = "/mnt/sda/gagelo01/Vcffile/Server_vcf/",
                                      vcffile_out,
                                      df_index,
                                      dt_gene_region,
                                      method_list = list("get_uni_cis", "get_coloc", "get_multicis")
                                     ) {

                                        
gwasvcf::set_bcftools()
gwasvcf::set_plink()

gene_region <- dt_gene_region[trait == df_index[id == ID, trait]]$gene_region

vcffile_exp <- paste0(out_wd, ID, "/", ID, ".vcf.gz")

res_all<- GagnonMR::run_all_pqtl_analyses(vcffile_exp = vcffile_exp, vcffile_out = vcffile_out, chrompos = gene_region, method_list = method_list)
setDT(res_all)
res_all[, study := df_index[id == ID, ]$consortium ]
return(res_all)
                                        
}

ID <- df_index[grepl("prot-2-|prot-3-", id)]$id %>% unique

run_all_pqtl_analyses_rosmap_banner_safely <- safely(run_all_pqtl_analyses_rosmap_banner)

options(future.globals.maxSize= 5e9)
plan(multicore, workers = 30, gc = TRUE) #I should try using multicore

tic()

res_all <-  future_map(as.list(ID)[1:40], function(x) {
  run_all_pqtl_analyses_rosmap_banner_safely(ID = x , #ensembl = ensembl_human,
                                             vcffile_out = "/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-1-1/trait-1-1.vcf.gz",
                                             method_list = list("get_uni_cis", "get_coloc", "get_multicis", "get_susie_coloc"),
                                            df_index = df_index_copy,
                                             dt_gene_region = dt_gene_region_copy)}, .options = furrr_options(seed = TRUE))

saveRDS(res_all, "Data/Modified/safely_res_all5.rds")

message("This script finished without errors")
toc()  
  


