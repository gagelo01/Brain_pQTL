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
somavec <- df_index[pmid == 34239129, trait]
dt_aptamers <- as.data.table(readat::aptamers)
all( somavec %in% dt_aptamers$AptamerId)
dtpath <- dt_aptamers[AptamerId %in% somavec, .(AptamerId, UniProt, EntrezGeneSymbol),]
vecpos <- GagnonMR::from_genecard_to_generegion(dtpath$EntrezGeneSymbol, window = 2e5)
vecpos[sapply(vecpos, function(x) length(x)>1)]<-NA #those that have two differents chromosome transform in NA
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

gene_region <- dt_gene_region[id == ID, ]$gene_region

vcffile_exp <- paste0(out_wd, ID, "/", ID, ".vcf.gz")

res_all<- GagnonMR::run_all_pqtl_analyses(vcffile_exp = vcffile_exp, vcffile_out = vcffile_out, chrompos = gene_region, method_list = method_list)
setDT(res_all)
res_all[, study := dt_gene_region[id == ID, study] ]
return(res_all)
                                        
}

ID <- df_index[grepl("prot-2-|prot-3-|prot-6-", id)]$id %>% unique
set.seed(seed = 2312)
ID <- sample(x = ID,size = length(ID),replace = FALSE)

run_all_pqtl_analyses_rosmap_banner_safely <- safely(run_all_pqtl_analyses_rosmap_banner)

options(future.globals.maxSize= 5e9)
plan(multicore, workers = 30, gc = TRUE) #plan(sequential)

top_cis_snp <- function( ID,
                          out_wd = "/mnt/sda/gagelo01/Vcffile/Server_vcf/",
                          dt_gene_region, 
                          ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"){
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  gene_region <- dt_gene_region[id == ID,]$gene_region
  vcffile_exp <- paste0(out_wd, ID, "/", ID, ".vcf.gz")
  dat_vcf <- gwasvcf::query_gwas(vcffile_exp, chrompos = gene_region)
  dat_tsmr <- gwasglue::gwasvcf_to_TwoSampleMR(dat_vcf)
  data.table::setDT(dat_tsmr)
  res <- dat_tsmr[which.min(pval.exposure), ]
  res <- res[,.(exposure, SNP, chr.exposure, pos.exposure, other_allele.exposure, effect_allele.exposure, eaf.exposure, 
         beta.exposure, se.exposure, pval.exposure)]
  res <- merge(dt_gene_region[id == ID, .(trait,UniProt,hgnc,study)], res, by.x = "trait", by.y = "exposure")
  res[,trait:=NULL]
  return(res) }
top_cis_snp_safely <-safely(top_cis_snp)
tic()

res_all <-  future_map(as.list(ID), function(x) {
  run_all_pqtl_analyses_rosmap_banner_safely(ID = x , #ensembl = ensembl_human,
                                             vcffile_out = "/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-1-1/trait-1-1.vcf.gz",
                                             method_list = list("get_uni_cis", "get_coloc", "get_multicis", "get_susie_coloc"),
                                            df_index = df_index_copy,
                                             dt_gene_region = dt_gene_region_copy)}, .options = furrr_options(seed = TRUE))

saveRDS(res_all, "Data/Modified/safely_res_all5.rds")

toc() 
top_cis <-  future_map(as.list(ID), function(x) {
  top_cis_snp_safely(ID = x , 
                     dt_gene_region = dt_gene_region_copy)}, .options = furrr_options(seed = TRUE))
saveRDS(top_cis, "Data/Modified/top_cis.rds")
fwrite(dt_gene_region, "Data/Modified/dt_gene_region.txt")
message("This script finished without errors")

  


