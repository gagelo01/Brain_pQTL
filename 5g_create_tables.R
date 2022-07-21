#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)
library(furrr)
library(gwasvcf)
library(TwoSampleMR)
library(xlsx)
library(writexl)

setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL")
res_pwas <- fread( "Data/Modified/res_pwas.txt")
res_map <- fread( "Data/Modified/res_map.txt")
data_united <- fread("Data/Modified/data_united_tpm.txt")
hypr_res <- fread( "Data/Modified/hypr_res.txt")
hypr_res <- hypr_res %>% distinct
ztest_ukb_giant <- fread( "Data/Modified/ztest_ukb_giant.txt")
susiesign <- fread("Data/Modified/susiesingstringent.txt" )
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")

# dt_res_hypr <- fread( "Data/Modified/dt_hypr_dietaryhabits.txt")
rescomparepqtl <- fread( "Data/Modified/ztest_rosmap_banner.txt")
# grs <- fread( "Data/Modified/grs_results_norm.txt")
# tfeq <- xlsx::read.xlsx("Data/Raw/TFEQ_Description.xlsx", sheetIndex = 1)
top_cis <- readRDS( "Data/Modified/top_cis.rds")
top_cis <- lapply(top_cis[sapply(top_cis, function(x) is.null(x$error))], function(x) x$result) %>% rbindlist(.,fill = TRUE)
top_cis <- top_cis[pval.exposure < 5e-8,]
#res_map
res_map_xlsx<- res_map[,.(UniProt, hgnc_symbol, study, posprob_coloc_PPH4, posprob_colocH4.SNP, posprob_colocH4.SNPexplained_var)]

#ztest 
ztest_ukb_giant <- merge(ztest_ukb_giant, dt_gene_region[,.(trait, hgnc,  study, UniProt)], by.x = c("exposure", "study") ,by.y = c("trait", "study"))
ztest_ukb_giant <- ztest_ukb_giant[, .(UniProt, hgnc, study, SNP,  z_stat, z_stat_pval, snpnotpresentingiant)]

colpwas<- c("UniProt", "hgnc_symbol", "study", colnames(res_pwas)[grepl("posprob_coloc|wald", colnames(res_pwas))])

colsusie<-colnames(susiesign)[grepl("multi_cis_susie|susiecoloc", colnames(susiesign))]
colsusie <- c("exposure",  "study", colsusie)

data_united_xlsx <-dcast(data_united, Description + fTau ~ tissue, value.var = "Averaged.TPM")
setnames(data_united_xlsx, "Description", "Gene")

setnames(hypr_res, paste0("trait", 1:3), c("trait1 (gene expression)", "trait2(protein levels)","trait3 (BMI)"))
hypr_res <- hypr_res[,.(gene_investigated, traits, posterior_prob, regional_prob, candidate_snp, posterior_explained_by_snp)]

data_united_xlsx <- data_united_xlsx[order(-fTau)]

# grs[, GRS := ifelse(GRS == "wGRS_Eloi_susie", "Brain_pQTL_GRS", "Top_BMI_GRS")]
# keycol <-c("GRS","p")
# grsa<- grs[!(DV %in% tfeq$TFEQ),]
# setorderv(grsa, keycol)
# keycol <-c("GRS","p")
# grsb<- grs[(DV %in% tfeq$TFEQ),]
# setorderv(grsb, keycol)

table_legend <- data.table(tables = c("ST1 Mapping results", "ST2 top p<5e-8 cis-snp sumstat",
                                      "ST3 uni-cis PWMR", 
                                      "ST4 ztest ukb and giant",
                                      "ST5 study heterogene.. lead snp",
                                      "ST6 multi-cis PWMR",
                                      "ST7 hyprcoloc RNA levels",
                                      "ST8 tpm per tissue"),
                           table_legend = c("Results of the mapping phase", "Lead cis SNPs for all genes and their relevant summary statistics. (only SNPs with p<5e-8 are shown)",
                                            "Result of the uni-cis protein-wide Mendelian randomization approach",
                                            "Ztest between the top cis SNPs of the PWMR approach for BMI from two different studies. The summary statistics from the Giant consortium and the summary statistics from the UKB",
                                            "Cochran's Q test between the lead SNPs of the ROS/MAP Banner and Yang protein datasets",
                                            "Result of the multi-cis protein-wide Mendelian randomization approach",
                                            "Results of the hyprcoloc analyses between RNA expression, protein expression and BMI",
                                            "Gene expression measure in transcript per million for all causal genes and for non-sexual tissues"),
                           column_legend = c("posprob_coloc_PPH4 = posterio probability of sharing the same causal varian; posprob_colocH4.SNP = The SNPs that is prioritised by coloc; posprob_colocH4.SNPexplained_var = The prioritised SNP explained variance of the colocalisation",
                                            "Column follows the format of the TowSampleMR R package",
                                            "b.wald, se.wald, pval.wald = statistics of the wald ratio; lead_snp.wald, rsq.exposure.wald, fstat.exposure.wald and pval_exposure.wald = statistics of the genetic instrument; posprob_coloc... = results of coloc",
                                            "snpnotpresentingiant = some snps were not present in giant (TRUE) therefore the z tests could not be calculated",
                                            "beta.exposure_study = the lead cis-SNP effect for the study; se.exposure_study = the lead cis-SNP standard error for the study; Q,p-value,df = Cochran's Q statistics.",
                                            "susiecoloc.n.credibleset.exposure susiecoloc.n.credibleset.outcome = the number of credible set prioritized by SuSiE; nsnpsusiecoloc = the number of SNP that colocalized at pph4>0.8",
                                            "columns follow the output in hyprcoloc R package",
                                            "each column is one tissue, each row is one gene"))

#
writexl::write_xlsx(x = list("Tables legend" = table_legend,
  "ST1 Mapping results" = res_map_xlsx,
                             "ST2 top p<5e-8 cis-snp sumstat" = top_cis,
                             "ST3 uni-cis PWMR" = res_pwas[, ..colpwas],
                             "ST4 ztest ukb and giant" = ztest_ukb_giant,
                             "ST5 study heterogene.. lead snp" = rescomparepqtl,
                             "ST6 multi-cis PWMR" = susiesign[,..colsusie],
                             "ST7 hyprcoloc RNA levels" = hypr_res,
                             "ST8 tpm per tissue" = data_united_xlsx
                             # "ST8-A GRS on dietary intake" = grsa,
                             # "ST8-B GRS on TFEQ" = grsb,
                             # "ST9 hyprcoloc dietary habits" = dt_res_hypr
                             ),
                    path = "Results/supplementary_tables.xlsx")

message("This script finished without errors")

