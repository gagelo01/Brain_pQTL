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
res_coloc <- readRDS( "Data/Modified/res_coloc_dh")
index <- sapply(res_coloc, function(x) is.null(x$result))
the_res <- lapply(res_coloc[!index], function(x) x$result) %>% rbindlist(., fill = TRUE)
the_res[, fdr.wald := p.adjust(pval.wald, method = "BH"), by = c("outcome", "study") ]
the_res <- the_res[posprob_coloc.mr > 0.8, ][fdr.wald < 0.05, ][steiger_dir.wald == TRUE][fstat.exposure.wald > 10,]

the_res[,split_index := paste0(exposure, "_", study)]
the_res_split <- split(the_res,  the_res$split_index)
saveRDS(the_res_split, "Data/Modified/the_res_split")

###
get_hypr_vcf <- function(list_vcffile, chrompos) {
  o <- gwasvcf::vcflist_overlaps(list_vcffile, chrompos = chrompos)
  veclength <- sapply(o, function(x) length(x))
  if(0 %in% veclength) {
    message("No overlaps in vcf1")
    return(NULL)
  }
  stopifnot(length(unique(veclength)) == 1)
  
  tab <- lapply(o, function(x) x %>%  gwasglue::gwasvcf_to_TwoSampleMR(., "exposure")) %>% rbindlist(., fill = TRUE)
  res <- GagnonMR::run_hypr_on_aligned(tab)
  return(res)
}

obtain_res_hypr <- function(res, df_index){
chrompos <- GagnonMR::from_genecard_to_generegion( res$exposure[1])
id_exp <- df_index[trait %in% res$exposure & consortium %in% res$study , ]$id
id_out <- df_index[pmid == 32193382 & trait %in% res$outcome,]$id
id_bmi <- "trait-1-1"
list_vcffile <- paste0("/mnt/sdf/gagelo01/Vcffile/Server_vcf/", c(id_exp, id_out, id_bmi), "/", c(id_exp, id_out, id_bmi), ".vcf.gz") %>% as.list
res_hypr <- get_hypr_vcf(list_vcffile = list_vcffile, chrompos = chrompos)
res_hypr[, gene_investigated :=  res$exposure[1]]
res_hypr[, study := res$study[1]]
return(res_hypr)
}


dt_res_hypr <- lapply(the_res_split, function(x) obtain_res_hypr(res = x, df_index = df_index)) %>% rbindlist(., fill = TRUE)


fwrite(dt_res_hypr, "Data/Modified/dt_hypr_dietaryhabits.txt")
