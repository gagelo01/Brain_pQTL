#!/usr/bin/env Rscript
library(GagnonMR)
library(tidyverse)
library(data.table)

setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL")

# load objects
res_all <- readRDS( "Data/Modified/safely_res_all5.rds")
index <- sapply(res_all, function(x) !is.null(x$result))
dt_res <- lapply(res_all[index], function(x) x$result) %>% rbindlist(., fill = TRUE)
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")

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
  dt_res[(chr.exposure.wald == region_df[i, ]$chr) & 
           (pos.exposure.wald >= region_df[i, ]$start) & (pos.exposure.wald <= 
                                                       region_df[i, ]$end), is_in_pleiotropic_region := TRUE]
}


#modify dt_res to include hgnc symbol
dt_aptamers <- as.data.table(readat::aptamers)
dt_res <- merge(dt_res,dt_aptamers[, .(AptamerId,UniProt) ], by.x = "exposure", by.y = "AptamerId", all.x = TRUE)
dt_res[study != "yang", UniProt := exposure]
dt_res <- merge(dt_res, dt_gene_region[, .(trait, id, gene_region, hgnc, ensg, study)],
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
res_map <- dt_res[posprob_coloc_PPH4 > 0.8 & hgnc_symbol %in% index_exposure_map, ] #results of the mapping

#results of the pwas
ntest <- dt_res[!is.na(b.wald), ]$UniProt %>% unique %>% length #number of different proteins tested number of test
uniprot <- dt_res[pval.wald < 0.05/ntest & posprob_coloc_PPH4 > 0.8  & steiger_pval.wald < 0.05 & fstat.exposure.wald > 10,]$UniProt
res_pwas <- dt_res[UniProt %in% uniprot,]

#####Replicatoin rate among uni-cis pwmr
k1 <- res_pwas[, .N, by = "UniProt"] #replicated
setnames(k1, "N", "present")
k2 <- res_pwas[pval.wald < 0.05/ntest & posprob_coloc_PPH4 > 0.8  & steiger_pval.wald < 0.05 & fstat.exposure.wald > 10, .N, by = "UniProt"] #replicated
setnames(k2, "N", "replicate")
k <- merge(k1,k2, by = "UniProt")
k[,ratio := replicate/present]
k[present == 1,category := "could not replicate"]
k[present>1 & ratio ==1, category := "replicated"]
k[present>1 & ratio<1, category := "not replicated"]
dt_replication_brainpqtl <- k

####top cis
gwasvcf::set_bcftools()
gwasvcf::set_plink()
arguments <- dt_res[!is.na(b.wald),.(id.exposure, lead_snp.wald)]
arguments <- arguments[lead_snp.wald!="",]
options(future.globals.maxSize= 5e9)
plan(multisession, workers = 20, gc = TRUE) #plan(sequential)

top_cis <-  future_map(split(arguments, 1:nrow(arguments)), function(x) {
  dat_tsmr <- gwasvcf::query_gwas(vcf = paste0("/mnt/sda/gagelo01/Vcffile/Server_vcf/", x$id.exposure, "/", x$id.exposure, ".vcf.gz"),
                                  rsid = x$lead_snp.wald) %>% gwasglue::gwasvcf_to_TwoSampleMR(.) %>% 
    as.data.table(.)
  return(dat_tsmr)}, .options = furrr_options(seed = TRUE)) %>% rbindlist(., fill=TRUE)

fwrite(top_cis,  "Data/Modified/top_cis.txt")
fwrite(dt_replication_brainpqtl, "Data/Modified/dt_replication_brainpqtl.txt")
fwrite(res_map, "Data/Modified/res_map.txt")
fwrite(res_pwas, "Data/Modified/res_pwas.txt")
fwrite(dt_res, "Data/Modified/dt_res_pwas.txt")

message("This script finished without errors")
