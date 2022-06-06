library(GagnonMR)
library(tidyverse)
library(data.table)
library(gwasvcf)
library(gwasglue)

# load objects
if(!file.exists("Data/Modified/dt_res.txt")) {
res_all <- readRDS( "Data/Modified/safely_res_all.rds")
index <- sapply(res_all, function(x) !is.null(x$result))
dt_res <- lapply(res_all[index], function(x) x$result) %>% rbindlist(., fill = TRUE)
} else {dt_res <- fread("Data/Modified/dt_res.txt")}
##########MAP###############

#select gene discovered by mapping
if(!file.exists("Data/Modified/BMI_loci.txt")) {
BMI_loci <- gwasvcf::query_gwas("/mnt/sdf/gagelo01/Vcffile/Server_vcf/trait-1-1/trait-1-1.vcf.gz", pval = 5e-8)
BMI_loci <- BMI_loci %>% gwasglue::gwasvcf_to_TwoSampleMR(., "outcome") %>% as.data.table
test <- ieugwasr::ld_clump(data.frame(rsid = BMI_loci$SNP,pval=BMI_loci$pval.outcome, id = BMI_loci$outcome),
                               clump_kb=10000,clump_r2=0.001,
                               plink_bin=genetics.binaRies::get_plink_binary(),
                               bfile="/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs")

BMI_loci <- BMI_loci[SNP %in% test$rsid,]
BMI_loci[,.N] #the number of independent genome-wide significant loci in yengo
fwrite(BMI_loci, "Data/Modified/BMI_loci.txt")
} else {
BMI_loci <- fread("Data/Modified/BMI_loci.txt")
}
# only include protein that are at +/- 100kb from BMI HITS
# obtain_generegion <- function(exp) {
# gene_region <- GagnonMR::from_genecard_to_generegion(exp, window = 0)
# 
# gene_region <- data.table::data.table( gene = exp,
# chr = sub(":.*", "", gene_region), 
# start = sub(".*:", "", gene_region) %>% sub("-.*", "", .),
# end = sub(".*:", "", gene_region) %>% sub(".*-", "", .))
# 
# return(gene_region)
# }
# 
# k <-  obtain_generegion(dt_res$exposure)
# k[, (c("chr", "start", "end")) := lapply(.SD, as.numeric), .SDcols = c("chr", "start", "end")]
# k <- distinct(k)
# list_exposure <- vector(mode = "list", length = BMI_loci[,.N])
# 
# for(i in 1:BMI_loci[,.N]) {
#   list_exposure[[i]]  <-  k[chr == BMI_loci[i, chr.outcome] & ((start <  BMI_loci[i, pos.outcome + 1e5] & start > BMI_loci[i, pos.outcome - 1e5])|
#                                                                  (end > BMI_loci[i, pos.outcome - 1e5] & end < BMI_loci[i, pos.outcome + 1e5])), unique(gene)]
# }
# 
# index <- sapply(list_exposure, length) > 1
# 
# test <- lapply(list_exposure[index], function(x) x %in% dt_res[posprob_coloc.mr > 0.8, unique(exposure)])
# any(sapply(test, function(x) sum(x)>1)) #It never returned more than one gene
# number <- sapply(test, function(x) length(x)>4 & sum(x)==1) %>% which
# list_exposure[index][number] 
# BMI_loci[index][number]
# index_exposure <- unlist(list_exposure) %>% unique
# 
# 
# 
# hist(unlist(list_loci))
# 
# 
# BMI_loci

dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")
gene_region <- dt_gene_region[, .(exposure, gene_region)] %>% distinct
gene_region <- separate(gene_region, col = "gene_region", into = c("chr", "start", "end"))
gene_region[, c("chr", "start", "end") := lapply(.SD, as.numeric), .SDcols = c("chr", "start", "end")]
gene_region_split <- split(gene_region, by = "exposure")
Nloci <- lapply(gene_region_split, function(bon) {
  Nloci <- BMI_loci[bon$chr == chr.outcome & bon$start < pos.outcome & pos.outcome < bon$end, .N]
  return(data.frame(exposure = bon$exposure, Nloci = Nloci))
}) %>% rbindlist

Nloci$Nloci %>% table

index_exposure <- Nloci[Nloci > 0,]$exposure %>% unique

res_map <- dt_res[posprob_coloc.mr > 0.8 & exposure %in% index_exposure, ] #results of the mapping
fwrite(res_map, "Data/Modified/res_map.txt")

###############results of the PWAS####################
#yang
res_yang <- fread( "Data/Modified/res_yang.txt")
res_yang[, study := "yang"]
res_yang[, exposure := gene.exposure]
res_yang[, gene.exposure := NULL]
dt_res <- rbindlist(list(res_yang, dt_res), fill = TRUE)

index_absent <- dt_res[, lapply(.SD, is.na), .SDcols = colnames(dt_res[, !c( "exposure",  "outcome", "study")])] %>% apply(., 1, function(x) all(x))
dt_res[index_absent, length(unique(exposure))] #for 0 proteins I could not map the gene to a gene region therefore I could not do analysis
dt_res <- dt_res[!index_absent, ]

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
ntest <- dt_res$exposure %>% unique %>% length #number of different proteins/different genes tested number of test
res_pwas <- dt_res[pval.wald < 0.05/ntest & (posprob_coloc.mr > 0.8 | ldcheck == TRUE) & pval_exposure.wald < 5e-8 & steiger_pval.wald < 0.05 & fstat.exposure.wald > 10,] 
fwrite(res_pwas, "Data/Modified/res_pwas.txt")
fwrite(dt_res, "Data/Modified/dt_res_pwas.txt")

###susiesign
colnames(dt_res)
nsnpcoloc <-dt_res[,  apply(.SD, 1, function(x) sum(x>0.8, na.rm = TRUE)), .SDcols = colnames(dt_res)[grepl("PP.H4.abf", colnames(dt_res))], by = c("exposure", "study")][order(exposure)]
setnames(nsnpcoloc, "V1", "nsnpsusiecoloc")
dt_res <- merge(dt_res, nsnpcoloc, by = c("exposure", "study"))
susiesign <- dt_res[pval.multi_cis_susie < 0.05/ntest & nsnpsusiecoloc > 1 & cochranQpval.multi_cis_susie > 0.05
                    & nsteigerfalse.multi_cis_susie == 0,]
susiesign[, .SD, .SDcols = colnames(susiesign)[grepl("susiecoloc.hit", colnames(susiesign))]]
susiesign <-susiesign[nsnpsusiecoloc == nsnp.multi_cis_susie,]

susiesign <- susiesign[minpvalexposure.multi_cis_susie < 1e-3,] #with at least one instrument with pval < 1e-3

#R2 LD between the genetic instruments
ldr2 <- apply(susiesign, 1 ,function(x) {
  ldmat <- ieugwasr::ld_matrix_local(c( x[names(x)=="susiecoloc.hit_1"],   x[names(x)=="susiecoloc.hit_2"]),
                                     plink_bin = genetics.binaRies::get_plink_binary(), 
                                     bfile = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs")
  return(data.frame(exposure = x[names(x)=="exposure"], study = x[names(x)=="study"], r2ld.multi_cis_susie = ldmat[1,2]^2))}) %>%
  rbindlist(., fill = TRUE)

susiesign<-merge(susiesign, ldr2, by=c("exposure", "study"))
fwrite(susiesign, file = "Data/Modified/susiesingstringent.txt" )
