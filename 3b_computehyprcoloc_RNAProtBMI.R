#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(gwasvcf)
library(gwasglue)
library(gwasglue)
library(GagnonMR)

setwd("/mnt/sde/gagelo01/Projects/Brain_pQTL")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
#load object
res_pwas <- fread( "Data/Modified/res_pwas.txt")
df_index <- fread("/mnt/sdf/gagelo01/Vcffile/server_gwas_id.txt")
res_map <- fread( "Data/Modified/res_map.txt")
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")

#define gene_to_include
arguments <- distinct(rbind(res_pwas[,.(exposure,study)], res_map[,.(exposure,study)]))
# gene_vec <- c(res_map$exposure, res_pwas$exposure) %>% unique
# arguments <- crossing(data.frame(exposure = gene_vec), study = c("rosmap", "banner"))
distinct(dt_gene_region[exposure %in% arguments$exposure, .(exposure, gene_region)])[,.N] == length(unique(arguments$exposure))

#brain eqtl
exposures_gtex <- fread("Data/Modified/exposures_gtex_hyprcoloc.txt")
#brain pqtl

extract_brain_pqtl <- function(gene_to_include, study, dt_gene_region) {
list_dat_tsmr <- vector(mode = "list", length = length(gene_to_include))
study_copy <- study

gene_region <- dt_gene_region[exposure == gene_to_include & study == study_copy,gene_region]
ID<- df_index[trait == gene_to_include & consortium == study, id]
vcffile_exp <- paste0("/mnt/sdf/gagelo01/Vcffile/Server_vcf/", ID, "/", ID, ".vcf.gz")
dat_tsmr <- query_gwas(vcffile_exp, chrompos = gene_region) %>% gwasvcf_to_TwoSampleMR(.) %>% as.data.table(.)
dat_tsmr[, exposure := paste0(study, "_", gene_to_include)]
dat_tsmr[, id.exposure := study]
dat_tsmr[, gene.exposure := gene_to_include]

return(dat_tsmr)
}

arguments_list <- split(arguments, seq(nrow(arguments)))
brain_pqtl_hypercoloc <- lapply(arguments_list, function(x)
  extract_brain_pqtl(gene_to_include = x$exposure, study = x$study, dt_gene_region = dt_gene_region )) %>%
  rbindlist(.,fill = TRUE)

#trait

list_dat_tsmr <- vector(mode = "list", length = length(gene_to_include))
for(i in 1:length(gene_to_include)) {
  gene_region <- dt_gene_region[exposure == gene_to_include[i], unique(gene_region)]
  ID<-"trait-1-1"
  vcffile_exp <- paste0("/mnt/sdf/gagelo01/Vcffile/Server_vcf/", ID, "/", ID, ".vcf.gz")
  dat_tsmr <- query_gwas(vcffile_exp, chrompos = gene_region) %>% gwasvcf_to_TwoSampleMR(.) %>% as.data.table(.)
  dat_tsmr[, id.exposure := exposure]
  dat_tsmr[, gene.exposure := gene_to_include[i]]
  dat_tsmr[, exposure := paste0(id.exposure, "-", gene.exposure)]
  list_dat_tsmr[[i]] <- dat_tsmr
}
all_exposures_yengo <- rbindlist(list_dat_tsmr)


list_hypercoloc <- list(eqtl = exposures_gtex, brain_pqtl = brain_pqtl_hypercoloc, trait = all_exposures_yengo)
dt_hyprcoloc <- rbindlist(list_hypercoloc, fill = TRUE)

fwrite(dt_hyprcoloc, "Data/Modified/dt_hyprcoloc.txt")


##Analyse
dt_hyprcoloc <- dt_hyprcoloc %>% distinct
all_id <- dt_hyprcoloc$id.exposure %>% unique 

id_to_include_list <- list(all_id[c(1,2,4)], all_id[c(1,3,4)])


obtain_list_align <- function(dt_hyprcoloc, id_to_include) {
list_hypr <- split(dt_hyprcoloc[(id.exposure %in% id_to_include),], by = "gene.exposure")
# list_hypr<-list_hypr[ sapply(list_hypr, function(x) length(unique(x$id.exposure)) == 3)]
list_aligned <- lapply(list_hypr, function(x) prepare_for_mvmr(x, x, harmonise_strictness = 1, should_clump = FALSE))
index <- sapply(list_aligned, function(x) x[,length(unique(id.exposure)) == length(id_to_include)])
list_aligned <- list_aligned[index]
return(list_aligned)
}


run_hypr_and_format <- function(list_align) {

  res <- map(list_align,  run_hypr_on_aligned)
  res_hypr <- rbindlist(res, fill = TRUE)
  
  res_hypr[, trait1 := ifelse(grepl("Brain_Frontal_Cortex_BA9", traits), TRUE, FALSE)]
  res_hypr[, trait2 := traits %>% ifelse(grepl("rosmap", .), TRUE, .) %>% 
             ifelse(grepl("banner", .), TRUE, .) ]
  res_hypr[trait2 != TRUE, trait2 := FALSE]
  # res_hypr[, trait3 :=ifelse(grepl("Whole_Blood", traits), TRUE, FALSE)]
  # res_hypr[, trait4 := ifelse( grepl("sun", traits), TRUE, FALSE)]
  res_hypr[,trait3 := ifelse(grepl("bmi_ukbgiant", traits), TRUE, FALSE)]
  res_hypr <- res_hypr[iteration == 1, ]
  
  bon <- res_hypr[, c("gene_investigated", paste0("trait", 1:3))]
  bon[, combination := paste(colnames(bon[, paste0("trait", 1:3)])[c(trait1==TRUE, trait2==TRUE, trait3==TRUE)], collapse = "-"), by = "gene_investigated"]
  mirge<- merge(res_hypr, bon[, .(gene_investigated, combination)], by = "gene_investigated")
  return(mirge)
}

list_align <- map(id_to_include_list, function(x) obtain_list_align(dt_hyprcoloc = dt_hyprcoloc, id_to_include = x))
list_align <- unlist(list_align, recursive = FALSE)
hypr_res <- run_hypr_and_format(list_align = list_align)
hypr_res[grepl("rosmap", allexposures), study := "rosmap"]
hypr_res[grepl("banner", allexposures), study := "banner"]
hypr_res <- hypr_res[order(gene_investigated, study)]
hypr_res[, paste0("trait",1:3) := lapply(.SD, as.logical), .SDcols = paste0("trait",1:3)]
hypr_res[, combination := apply(.SD, 1, function(x) paste(paste0("trait",1:3)[c(x[1]==TRUE, x[2]==TRUE, x[3]==TRUE)], collapse = "-") ), .SDcols = paste0("trait",1:3)]
hypr_res <- distinct(hypr_res)
hypr_res[, exposure_study := paste0(gene_investigated, "_", study)]
fwrite(hypr_res, "Data/Modified/hypr_res.txt")
saveRDS(list_align, "Data/Modified/list_align")
message(paste0("script finised without errors at ", Sys.time()))
