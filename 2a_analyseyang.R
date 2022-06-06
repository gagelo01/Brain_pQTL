#!/usr/bin/env Rscript
library(data.table)
library(TwoSampleMR)
library(tidyverse)
library(gwasvcf)
library(gwasglue)
library(GagnonMR)

setwd("/mnt/sde/gagelo01/Projects/Brain_pQTL")
gwasvcf::set_bcftools()
gwasvcf::set_plink()

inst_yang <- fread( "Data/Modified/inst_yang_clump.txt")
bmi_out <- gwasvcf::query_gwas(vcf = "/mnt/sdf/gagelo01/Vcffile/Server_vcf/trait-1-1/trait-1-1.vcf.gz", rsid = inst_yang$SNP) %>%
  gwasglue::gwasvcf_to_TwoSampleMR(., type = "outcome")


ldcheck_yang <- function(exp) {
  message(paste0("**********Initialising  ", exp, "**************"))
  dt_exposure <- inst_yang[exposure == exp,]
  window <- (2e6)/2
  snp_region <- dt_exposure[, paste0( chr.exposure, ":", (pos.exposure-window) %>% ifelse(.<1,"1",.), "-", pos.exposure+window)]
  all_out <- gwasvcf::query_gwas("/mnt/sdf/gagelo01/Vcffile/Server_vcf/trait-1-1/trait-1-1.vcf.gz", chrompos = snp_region)
  dt_outcome <- all_out %>% gwasglue::gwasvcf_to_TwoSampleMR(., type = "outcome") %>% as.data.table(.)
  dt_outcome30 <- dt_outcome[order(pval.outcome),][1:30]
  
  snp_proxy_vec <- tryCatch(
    expr = {gwasvcf::get_ld_proxies(rsid = dt_exposure$SNP,
                                    bfile = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs", tag_r2 = 0.8)$SNP_B    },
    error = function(e){return(list()) })
  
  df_res <-data.frame(exposure = exp, ldcheck = if(length(snp_proxy_vec) == 0) {NA}else{ dt_outcome30[, any(SNP %in% snp_proxy_vec)] })
  return(df_res)
}

res_ldcheck <- lapply(as.list(unique(inst_yang$exposure)), function(dat) ldcheck_yang(dat)) %>% rbindlist(., fill = TRUE)

#MR
harm <- TwoSampleMR::harmonise_data(inst_yang, bmi_out, 1)
harm <- TwoSampleMR::steiger_filtering( harm )
harm <- TwoSampleMR::add_rsq(harm)
harm$fstat.exposure <- fstat_fromdat(harm)
setDT(harm)
res_all <- TwoSampleMR::mr(harm, method_list = c("mr_wald_ratio", "mr_ivw")) %>% data.table::as.data.table(.)
res_all[, c("id.exposure", "id.outcome", "method") :=NULL]
res_wald <- merge(res_all[nsnp ==1,], harm[,.(exposure, SNP, chr.exposure, pos.exposure, pval.exposure, steiger_dir, steiger_pval, rsq.exposure, fstat.exposure )], by = "exposure")

old = c("b","se", "pval", "SNP", "chr.exposure", "pos.exposure", "pval.exposure", "steiger_dir", "steiger_pval", "rsq.exposure", "fstat.exposure")
new = c(paste0(c("b", "se", "pval"), ".wald"), "lead_snp.wald", "chr.exposure.wald", "pos.exposure.wald", "pval_exposure.wald", "steiger_dir.wald",
"steiger_pval.wald",  "rsq.exposure.wald", "fstat.exposure.wald")
data.table::setnames(res_wald, old = old, new = new)

res_ivw <- res_all[nsnp >1,]
setnames(res_ivw, c("b","se", "pval"), c("b.ivw", "se.ivw", "pval.ivw"))
res_ivw[,nsnp := NULL]
res_wald <- rbindlist(list(res_wald, res_ivw), fill = TRUE)
#merge
res_yang <- merge(res_wald, res_ldcheck, by = c("exposure"))
res_yang <- merge(res_yang, inst_yang[, .(exposure, gene.exposure)], by = "exposure")
fwrite(res_yang, "Data/Modified/res_yang.txt")
message("this script finished without errors")
