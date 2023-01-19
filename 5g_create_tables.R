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
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")
dt_replication_brainpqtl <- fread( "Data/Modified/dt_replication_brainpqtl.txt")
coloc_matrix <- readRDS( "Data/Modified/coloc_matrix.rds")
dt_res <- fread( "Data/Modified/dt_res_pwas.txt")
res_specificity <- fread("Data/Modified/res_specificity.txt")
top_cis <- fread("Data/Modified/top_cis.txt")

#res_map
res_map_xlsx<- res_map[,.(UniProt, hgnc_symbol, study, posprob_coloc_PPH4, posprob_colocH4.SNP, posprob_colocH4.SNPexplained_var)]

#top_cis
top_cis[, c("ncase.exposure", "ncontrol.exposure", "id.exposure"):=NULL]
#ztest 
ztest_ukb_giant <- merge(ztest_ukb_giant, dt_gene_region[,.(trait, hgnc,  study)], by.x = c("exposure", "study") ,by.y = c("trait", "study"))
ztest_ukb_giant <- ztest_ukb_giant[, .(exposure, UniProt, hgnc, study, SNP,  z_stat, z_stat_pval, snpnotpresentingiant)]

colpwas<- c("UniProt", "hgnc_symbol", "study", colnames(res_pwas)[grepl("posprob_coloc|wald", colnames(res_pwas))])

#matrix coloc
coloc_proportion <- as.data.frame(coloc_matrix$mat_proportion)
coloc_proportion<- cbind(" "=rownames(coloc_proportion), coloc_proportion)
setDT(coloc_proportion)
#res_specificity
setcolorder(res_specificity, c())

#hypr_res
hypr_res <- hypr_res[, .SD[1], by = "gene_investigated"]

#data tpm
data_united_xlsx <-dcast(data_united, Description + fTau ~ tissue, value.var = "Averaged.TPM")
setnames(data_united_xlsx, "Description", "Gene")

hypr_res <- hypr_res[,.(gene_investigated, traits, posterior_prob, regional_prob, candidate_snp, posterior_explained_by_snp)]

data_united_xlsx <- data_united_xlsx[order(-fTau)]

######cleanify
list_all_tables <- list("ST1 Mapping results" = res_map_xlsx,
     "ST2 top p<5e-8 cis-snp sumstat" = top_cis,
     "ST3 uni-cis PWMR" = res_pwas[, ..colpwas],
     "ST4 ztest ukb and giant" = ztest_ukb_giant,
     "ST5 replication rate" = dt_replication_brainpqtl,
     "ST6 replication other tissues" = res_specificity,
     "ST7 colocalisation matrix ratio" = coloc_proportion,
     "ST8 hyprcoloc RNA levels" = hypr_res,
     "ST9 tpm per tissue" = data_united_xlsx
)

format_table_number <- function(x) {
  if(is.numeric(x)){
    if(any(x[!is.na(x)]>=10000 | x[!is.na(x)]<0.001)) {
      formatC(x, format = "e", digits = 1)
    } else(round(x = x, digits = 3))
  } else {x}
}

list_all_tables <- map(list_all_tables, function(x) x[, lapply(.SD, format_table_number)])

table_legend <- data.table(tables = c("ST1 Mapping results",
                                      "ST2 top p<5e-8 cis-snp sumstat",
                                      "ST3 uni-cis PWMR", 
                                      "ST4 ztest ukb and giant",
                                      "ST5 replication rate",
                                      "ST6 tpm per tissue",
                                      "ST7 replication other tissues",
                                      "ST8 colocalisation matrix ratio",
                                      "ST9 hyprcoloc RNA levels"),
                           table_legend = c("Results of the mapping phase", 
                                            "Lead cis SNPs for all genes and their relevant summary statistics. (only SNPs with p<5e-8 are shown)",
                                            "Result of the uni-cis protein-wide Mendelian randomization approach",
                                            "Ztest between the top cis SNPs of the PWMR approach for BMI from two different studies. The summary statistics from the Giant consortium and the summary statistics from the UKB",
                                            "Replication rate among the three brain pQTL datasets used",
                                            "Gene expression measure in transcript per million for all causal genes and for non-sexual tissues",
                                            "Result of the uni-cis protein-wide Mendelian randomization approach for other tissues",
                                            "Matrix : what is the percentage of proteins with pph4 > 0.8 accross studies",
                                            "Results of the hyprcoloc analyses between RNA expression, protein expression and BMI"),
                           column_legend = c("posprob_coloc_PPH4 = posterio probability of sharing the same causal varian; posprob_colocH4.SNP = The SNPs that is prioritised by coloc; posprob_colocH4.SNPexplained_var = The prioritised SNP explained variance of the colocalisation",
                                            "Column follows the format of the TowSampleMR R package",
                                            "b.wald, se.wald, pval.wald = statistics of the wald ratio; lead_snp.wald, rsq.exposure.wald, fstat.exposure.wald and pval_exposure.wald = statistics of the genetic instrument; posprob_coloc... = results of coloc",
                                            "snpnotpresentingiant = some snps were not present in giant (TRUE) therefore the z tests could not be calculated",
                                            "UniProt = the UniProt id, present = the number of dataset where the UniProt has been assessed, the number of study that confirm the protein is causal, ratio = ratio of column 2 and 3, category = a summary in three categories",
                                            "each column is one tissue, each row is one gene",
                                             "b.wald, se.wald, pval.wald = statistics of the wald ratio; lead_snp.wald, rsq.exposure.wald, fstat.exposure.wald and pval_exposure.wald = statistics of the genetic instrument; posprob_coloc... = results of coloc",
                                             "this is a matrix, where each cell represents the proportion of protein with a PPH$ > 0.8", 
                                           "columns follow the output in hyprcoloc R package: gene_investigated = the gene that was investigated. traits = the traits that were deemed to hyprcolocalise. "))

#

writexl::write_xlsx(x = c("Tables legend" = list(table_legend), list_all_tables),
                    path = "Results/supplementary_tables.xlsx")

message("This script finished without errors")

