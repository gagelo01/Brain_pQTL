#!/usr/bin/env Rscript
library(GagnonMR)
library(tidyverse)
library(data.table)

setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL")

# load objects
res_all <- readRDS( "Data/Modified/safely_res_all5.rds")
index <- sapply(res_all, function(x) !is.null(x$result))
dt_res <- lapply(res_all[index], function(x) x$result) %>% rbindlist(., fill = TRUE)

##########MAP###############

#select gene discovered by mapping
if(!file.exists("Data/Modified/BMI_loci.txt")) {
  BMI_loci <- gwasvcf::query_gwas("/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-1-1/trait-1-1.vcf.gz", pval = 5e-8)
  BMI_loci <- BMI_loci %>% gwasglue::gwasvcf_to_TwoSampleMR(., "outcome") %>% as.data.table
  test <- ieugwasr::ld_clump(data.frame(rsid = BMI_loci$SNP,pval=BMI_loci$pval.outcome, id = BMI_loci$outcome),
                             clump_kb=10000,clump_r2=0.001,
                             plink_bin=genetics.binaRies::get_plink_binary(),
                             bfile="/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs")
  
  BMI_loci <- BMI_loci[SNP %in% test$rsid,]
  BMI_loci[,.N] #the number of independent genome-wide significant loci in pulit
  fwrite(BMI_loci, "Data/Modified/BMI_loci.txt")
} else {
  BMI_loci <- fread("Data/Modified/BMI_loci.txt")
}
# only include protein that are at +/- 100kb from BMI HITS


###############results of the PWAS####################
#report if results are in the ABO, APOE, HLA-A gene region
setnames(dt_res, c("chr.exposure.wald","pos.exposure.wald"), c("chr.exposure", "pos.exposure"))
to_exclude = c("APOE", "ABO", "HLA-A")
window = 2e+06
gencode <- fread("/home/couchr02/Mendel_Commun/Nicolas/GTEx/gencode.v19.genes.v7.patched_contigs.txt")
list <- vector(mode = "list", length = length(to_exclude))
for (i in 1:length(to_exclude)) {
  bon <- gencode[gene_name == to_exclude[i], ]
  list[[i]] <- data.frame(chr = bon[1, ]$chr, start = min(bon$start) - 
                            window/2, end = max(bon$end) + window/2, gene_name = bon[1, 
                            ]$gene_name)
}
region_df <- rbindlist(list)
dt_res[, is_in_pleiotropic_region := FALSE]
for (i in 1:nrow(region_df)) {
  dt_res[(chr.exposure == region_df[i, ]$chr) & 
           (pos.exposure >= region_df[i, ]$start) & (pos.exposure <= 
                                                       region_df[i, ]$end), is_in_pleiotropic_region := TRUE]
}


#k
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")
dt_aptamers <- as.data.table(readat::aptamers)
dt_res <- merge(dt_res,dt_aptamers[, .(AptamerId,UniProt) ], by.x = "exposure", by.y = "AptamerId", all.x = TRUE)
dt_res[study != "yang", UniProt := exposure]
dt_res <- merge(dt_res, dt_gene_region[, .(trait, id, gene_region, hgnc, study)],
                by.x = c("exposure", "study"), by.y = c("trait", "study"), all.x = TRUE)
setnames(dt_res, "hgnc","hgnc_symbol")
##results of the mapping
gene_region <- separate(dt_gene_region, col = "gene_region", into = c("chr", "start", "end"))
gene_region[, c("chr", "start", "end") := lapply(.SD, as.numeric), .SDcols = c("chr", "start", "end")]
gene_region_split <- split(gene_region[!is.na(start),], by = "id")
Nloci <- lapply(gene_region_split, function(bon) {
  Nloci <- BMI_loci[bon$chr == chr.outcome & bon$start < pos.outcome & pos.outcome < bon$end, .N]
  return(data.frame(hgnc = bon$hgnc, Nloci = Nloci))
}) %>% rbindlist

Nloci$Nloci %>% table

index_exposure_map <- Nloci[Nloci > 0,]$hgnc %>% unique
res_map <- dt_res[posprob_coloc_PPH4 > 0.8 & UniProt %in% index_exposure_map, ] #results of the mapping

#results of the pwas
ntest <- dt_res$UniProt %>% unique %>% length #number of different proteins tested number of test
res_pwas <- dt_res[pval.wald < 0.05/ntest & posprob_coloc_PPH4 > 0.8 & pval_exposure.wald < 5e-8 & steiger_pval.wald < 0.05 & fstat.exposure.wald > 10,] 

###susiesign
nsnpcoloc <-dt_res[,  apply(.SD, 1, function(x) sum(x>0.8, na.rm = TRUE)), .SDcols = colnames(dt_res)[grepl("PP.H4.abf", colnames(dt_res))], by = c("exposure", "study")][order(exposure)]
setnames(nsnpcoloc, "V1", "nsnpsusiecoloc")
dt_res <- merge(dt_res, nsnpcoloc, by = c("exposure", "study"))

#multi-cis
susiesign <- dt_res[pval.multi_cis < 0.05/ntest & nsnpsusiecoloc > 0   & cochranQpval.multi_cis > 0.05,]
fwrite(res_map, "Data/Modified/res_map.txt")
fwrite(susiesign, file = "Data/Modified/susiesingstringent.txt" )
fwrite(res_pwas, "Data/Modified/res_pwas.txt")
fwrite(dt_res, "Data/Modified/dt_res_pwas.txt")

message("This script finished without errors")
