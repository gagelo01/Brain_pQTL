#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(GagnonMR)

gwasvcf::set_bcftools()
gwasvcf::set_plink()
ldref<-"/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs"
res_pwas <- fread( "Data/Modified/res_pwas.txt")
res_map <- fread( "Data/Modified/res_map.txt")
df_index <- fread("/mnt/sdf/gagelo01/Vcffile/server_gwas_id.txt")
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

ao <- fread("/mnt/sdf/gagelo01/Vcffile/available_outcomes_2021-10-13.txt")
ao[id %in% list.files("/mnt/sdf/gagelo01/Vcffile/MRBase_vcf/")]

idbmi <- c("ieu-a-835", "ukb-b-19953")
sumstat<- lapply(as.list(idbmi), function(x)
  gwasvcf::query_gwas(paste0("/mnt/sdf/gagelo01/Vcffile/MRBase_vcf/", x, "/", x, ".vcf.gz"), rsid = snp_z, proxies = "no") %>%
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
res_combine <- merge(res_combine[,.(exposure,  study, lead_snp.wald)], z_res, by.x = "lead_snp.wald", by.y = "SNP") %>% distinct
setnames(res_combine, "lead_snp.wald", "SNP")
fwrite(res_combine, "Data/Modified/ztest_ukb_giant.txt")



####banner rosmap
res_pwas[,c(setdiff(colnames(res_pwas), colnames(res_map))) := NULL]
res_map[, setdiff(colnames(res_map), colnames(res_pwas)) := NULL]
res_combine <- rbindlist(list(res_pwas, res_map), fill = TRUE)
res_combine <- distinct(res_combine)
args <- res_combine[, .(lead_snp.wald, exposure) ] %>% distinct
args_list <- split(args, 1:nrow(args))


comparestudyz <- function(exposure, SNP_z) {
ID <- df_index[trait == exposure,]$id
sumstat<- lapply(as.list(ID), function(x){
  k <- tryCatch(expr = {
    gwasvcf::query_gwas(paste0("/mnt/sdf/gagelo01/Vcffile/Server_vcf/", x, "/", x, ".vcf.gz"), rsid = SNP_z, proxies = "no", bfile = ldref) %>%
    gwasglue::gwasvcf_to_TwoSampleMR(.) %>%
    as.data.table  }, error = function(e) {
      return(data.table(exposure = exposure, SNP = SNP_z))
    })
  k[, study := df_index[id == x, consortium]]
  return(k)}) %>% rbindlist(.,fill=TRUE)

zres_wide<- dcast(sumstat, exposure ~ study, value.var = c("beta.exposure", "se.exposure"))

if(nrow(sumstat) == 1){return(cbind(zres_wide, presentinbothstudy = FALSE))}

z_res <-Z_test(mean.x = sumstat[1,beta.exposure], mean.y = sumstat[2,beta.exposure],se.x = sumstat[1,se.exposure], se.y = sumstat[2,se.exposure])
return(cbind(zres_wide, SNP= SNP_z, z_res, presentinbothstudy = TRUE))
}

comparestudyz_safely <- safely(comparestudyz)
rescompare <- map(args_list, function(x) comparestudyz(exposure = x$exposure, SNP_z = x$lead_snp.wald)) %>%
  rbindlist(., fill = TRUE) %>% as.data.table

rescompare[,z_stat_pval := NULL]
rescompare[,heterogenous_zstatabove5 := ifelse(abs(z_stat)>5, TRUE, FALSE)]
rescompare[heterogenous_zstatabove5==TRUE,]

fwrite(rescompare, "Data/Modified/ztest_rosmap_banner.txt")
