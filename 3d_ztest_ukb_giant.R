#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(GagnonMR)
setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL/")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref<-"/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
res_pwas <- fread( "Data/Modified/res_pwas.txt")
res_map <- fread( "Data/Modified/res_map.txt")
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")

# susiesign <- fread("Data/Modified/susiesingstringent.txt" )
res_pwas[,c(setdiff(colnames(res_pwas), colnames(res_map))) := NULL]
res_map[, setdiff(colnames(res_map), colnames(res_pwas)) := NULL]
res_combine <- rbindlist(list(res_pwas, res_map), fill = TRUE)
res_combine <- distinct(res_combine)
snp_z <- res_combine$lead_snp.wald %>% unique

Z_test<- function(mean.x, mean.y, se.x, se.y, Ho = 0) {
  stderr <- sqrt(se.x^2 + se.y^2)
  zstat <- (mean.x - mean.y - Ho)/stderr
  pval <- 2 * pnorm(-abs(zstat))
  return(data.frame(z_stat = zstat,
                    z_stat_pval = pval))
}


idbmi <- c("ieu-a-835", "ukb-b-19953")
sumstat<- lapply(as.list(idbmi), function(x)
  gwasvcf::query_gwas(paste0("/mnt/sda/gagelo01/Vcffile/MRBase_vcf/", x, "/", x, ".vcf.gz"), rsid = snp_z, proxies = "no") %>%
    gwasglue::gwasvcf_to_TwoSampleMR(.)) %>% rbindlist(.,fill=TRUE)

sumstat_split <- split(sumstat, sumstat$SNP)
z_res<- lapply(sumstat_split, function(x) {
zres <- Z_test(mean.x = x[exposure == "UKB-b-19953", beta.exposure], mean.y = x[exposure == "ieu-a-835", beta.exposure],
       se.x = x[exposure == "UKB-b-19953", se.exposure], se.y = x[exposure == "ieu-a-835", se.exposure])
if(nrow(zres) == 0) {zres<-data.frame(z_stat = NA, z_stat_pval = NA)}
return(cbind(SNP = x$SNP[1],zres))}) %>% rbindlist(.,fill = TRUE)



k<-sumstat[, paste(unique(exposure), collapse = "_"), by = "SNP"]
k[, snpnotpresentingiant := !grepl("ieu-a-835", V1)]
z_res <- merge(z_res,k[, .(SNP, snpnotpresentingiant)], by = "SNP" )
res_combine <- merge(res_combine[,.(exposure, UniProt,  study, lead_snp.wald)], z_res, by.x = "lead_snp.wald", by.y = "SNP") %>% distinct
setnames(res_combine, "lead_snp.wald", "SNP")


sumstat[, exposure:= exposure %>% gsub("UKB-b-19953", "UKB", .) %>% gsub("ieu-a-835", "GIANT", .)]
k <- sumstat[,.(exposure, SNP, beta.exposure, se.exposure)]
k<- data.table::dcast(k, SNP ~ exposure, value.var = c("beta.exposure", "se.exposure"))
res_combine <- merge(k, res_combine, by = "SNP")

fwrite(res_combine, "Data/Modified/ztest_ukb_giant.txt")

####banner rosmap
# args <- res_pwas[, .(lead_snp.wald, UniProt) ] %>% distinct
# args_list <- split(args, 1:nrow(args))
# 
# comparestudyz <- function(prot, SNP_z) {
# ID <- dt_gene_region[UniProt == prot,]$id
# sumstat<- lapply(as.list(ID), function(x){
#   k <- tryCatch(expr = {
#     gwasvcf::query_gwas(paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", x, "/", x, ".vcf.gz"), rsid = SNP_z, proxies = "no", bfile = ldref) %>%
#     gwasglue::gwasvcf_to_TwoSampleMR(.) %>%
#     as.data.table  }, error = function(e) {
#       return(data.table(id.exposure = x, SNP = SNP_z))
#     })
#   k[, study := df_index[id == x, consortium]]
#   k[, study := study %>% gsub("None", "Yang", .)]
#   k[,id.exposure := x]
#   return(k)}) %>% rbindlist(.,fill=TRUE)
# 
# sumstat <- merge(sumstat, dt_gene_region[,.(id,UniProt)], by.x = c("id.exposure"), by.y = c("id"))
# zres_wide<- dcast(sumstat, UniProt ~ study, value.var = c("beta.exposure", "se.exposure"))
# 
# if(nrow(sumstat) == 1){return(cbind(zres_wide, SNP = SNP_z, presentin_nstudy = 1))}
# 
# res <- mada::cochran.Q(x = sumstat$beta.exposure, weights = sumstat$se.exposure/1)
# 
# return(cbind(zres_wide, SNP= SNP_z, as.data.table(as.list(res)), presentin_nstudy = sumstat[,.N]))
# }
# 
# comparestudyz_safely <- safely(comparestudyz)
# rescompare <- map(args_list, function(x) comparestudyz(prot = x$UniProt, SNP_z = x$lead_snp.wald)) %>%
#   rbindlist(., fill = TRUE) %>% 
#   as.data.table(.)
# 
# fwrite(rescompare, "Data/Modified/ztest_rosmap_banner.txt")
message("This script finished without errors")