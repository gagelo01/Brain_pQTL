#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(GagnonMR)
gwasvcf::set_bcftools()
gwasvcf::set_plink()
setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL")
ldref <- "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"

res_map <- fread( "Data/Modified/res_map.txt")
res_pwas <- fread( "Data/Modified/res_pwas.txt")
susiesign <- fread("Data/Modified/susiesingstringent.txt" )
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")

#resuni
res_map[, SNP := posprob_colocH4.SNP]
res_pwas[, SNP := lead_snp.wald]
resuni <- rbind(res_map[,.(id,SNP,hgnc_symbol, UniProt,study)], res_pwas[,.(id,SNP,hgnc_symbol, UniProt, study)]) %>% distinct(.)
resuni[,approach := "uni"]
##the object to imput FUMA
resmulti <- melt(susiesign, id.vars = c("id","hgnc_symbol", "UniProt", "study"), measure.vars = colnames(susiesign)[grepl("susiecoloc.hit_", colnames(susiesign))])
resmulti <- resmulti[value != "",]
resmulti[,variable:=NULL]
resmulti[,approach := "multi"]
setnames(resmulti, "value", "SNP")
resall <- rbind(resmulti, resuni)
resall_split <- split(resall, resall$id)
res <- map(resall_split, function(x){
tsmr<- gwasvcf::query_gwas(vcf = paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", x$id, "/", x$id, ".vcf.gz"),
                        rsid = x$SNP, proxies = "no") %>% 
  gwasglue::gwasvcf_to_TwoSampleMR(.) %>% as.data.table(.)
tsmr <- merge(tsmr, x , by = "SNP")
return(tsmr)}) %>% rbindlist(.,fill = TRUE)

fwrite(res, "Data/Modified/input_fuma.txt")
message("This script finished without errors")

