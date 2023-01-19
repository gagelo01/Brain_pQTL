#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)
library(furrr)

setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL")
gwasvcf::set_bcftools()
gwasvcf::set_plink()

dt_tra <- fread("Data/Modified/dt_specificity_gene_region.txt")
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")
dt_tra <- rbind(dt_tra, dt_gene_region, fill = TRUE)
res_pwas <- fread( "Data/Modified/res_pwas.txt")
res_map <- fread( "Data/Modified/res_map.txt")
dt_res <- fread("Data/Modified/dt_res_pwas.txt")

######
parameters <- GagnonMR::default_param()
parameters$path <- c(default_param()$path, "/mnt/sda/gagelo01/Projects/Brain_pQTL/Data/Modified/Gtex_vcf/")
#######
message("Initialising hyprcoloc")
arguments <- dt_tra
arguments[,out_wd := ifelse(study == "GTEX", 
                            "/mnt/sda/gagelo01/Projects/Brain_pQTL/Data/Modified/Gtex_vcf/",
                            "/mnt/sda/gagelo01/Vcffile/Server_vcf/")]
arguments[,  vcfpath :=paste0(out_wd, id, "/", id, ".vcf.gz")]
arguments <- arguments[,.(id, vcfpath, UniProt, gene_region, hgnc, study)]

uniprottoloselect <- dt_tra[UniProt%in%c(res_pwas$UniProt,res_map$UniProt), .N, by = "UniProt"][N>2]$UniProt
arguments_hypr <- arguments[UniProt %in% uniprottoloselect,]
arguments_hypr <- arguments_hypr[study %in% c("GTEX", "yang", "rosmap",  "banner")] #Since this analysis only assess wheteher the same genetic variants also affect gene expression we remove data from the blood
arguments_hypr <- split(arguments_hypr, by = "hgnc")
options(future.globals.maxSize = 5e9)
plan(multisession, workers =10, gc = TRUE)

exposure_hypr <- future_map(arguments_hypr, function(x) { 
  gwasvcf::set_bcftools()
       res<- get_hyprcoloc(vcffile_exp = x$vcfpath, 
                     vcffile_out = "/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-1-1/trait-1-1.vcf.gz",
                     chrompos = unique(x$gene_region),
                     parameters = parameters)$align
       res$gene.exposure <- unique(x$hgnc)
       return(res)}, .options = furrr_options(seed = TRUE)) %>% 
  rbindlist(., fill = TRUE)

exposure_hypr <- merge(exposure_hypr, dt_tra[,.(id,study)], by.x = "id.exposure", by.y = "id", all.x = TRUE)
exposure_hypr[study%in% dt_tra$study,exposure:=paste0(study, "_", gene.exposure, "_", id.exposure)]

res_hypr <- map(split(exposure_hypr, exposure_hypr$gene.exposure), GagnonMR::run_hypr_on_aligned) %>% 
  rbindlist(., fill = TRUE)

########The question is among the evaluated proteins (proteins with pQTL) what is the colocalisation rate ?
message("Initialising coloc")

plan(multisession, workers =20, gc = TRUE)

uniprot_vec <- dt_res[!is.na(pval.exposure.wald), unique(UniProt)]
k<-arguments[, .N, by = "UniProt"][N>1]$UniProt
uniprot_vec <- uniprot_vec[uniprot_vec%in%k]
arguments_coloc <- split(arguments[UniProt %in% uniprot_vec,], by = "hgnc")

res_coloc <- future_map(arguments_coloc, function(x) {
  gwasvcf::set_bcftools()
  gwasvcf::set_plink()
  k <- combn(x$vcfpath, 2)
  list_vcf <- plyr::alply(k,2)
  res_coloc<- map(list_vcf, function(y) {
    GagnonMR::get_coloc(vcffile_exp = y[1], vcffile_out = y[2], chrompos = x[1,]$gene_region, parameters = parameters)
  }) %>% rbindlist(., fill = TRUE)
  return(res_coloc)}, .options = furrr_options(seed = TRUE)) %>%
  rbindlist(., fill = TRUE)

res_coloc <- merge(res_coloc, dt_tra[,.(id, UniProt, study, gene_region, hgnc)], by.x = "id.exposure", by.y = "id")
setnames(res_coloc, "study", "study_exposure")
res_coloc <- merge(res_coloc, dt_tra[,.(id, study)], by.x = "id.outcome", by.y = "id")
setnames(res_coloc, "study", "study_outcome")


res_coloc[, coloc := posprob_coloc_PPH4>0.8]
res_coloc[!(UniProt%in%res_pwas$UniProt),]
res_coloc[ study_exposure == "banner" & study_outcome == "rosmap", c("study_exposure", "study_outcome") := .(study_outcome, study_exposure)]
####create posprob matrix
message("Initialising coloc matrix")

format_test_to_matrix<- function(test) {
  setnames(test, c("study_exposure", "study_outcome", "V1"), c("exposure", "outcome","value"))
  test2 <- test
  colnames(test2)<-c("outcome", "exposure", "value")
  k<-c(test$exposure,test$outcome) %>% unique
  k <- k[!(k%in%test[exposure==outcome, exposure])]
  test3<-data.frame(exposure = k, outcome = k, value = NA)
  test4 <- rbindlist(list(test, test2, test3), use.names = TRUE, fill = TRUE)
  test4<-distinct(test4)
  setDT(test4)
  test4<-test4[order(exposure, outcome),]
  test4 <- data.table::dcast(test4, outcome ~ exposure, value.var = "value")
  test4 <-as.matrix(test4, rownames = "outcome")
}

test <- res_coloc[ , sum(coloc, na.rm = TRUE)/.N, by = c("study_exposure", "study_outcome")]
mat_ratio <- format_test_to_matrix(test)
otter_dendro <- as.dendrogram(hclust(d = dist(x = mat_ratio)))
otter_order <- order.dendrogram(otter_dendro)
otter_order<-rownames(mat_ratio)[otter_order]
mat_ratio<-mat_ratio[otter_order, otter_order]

test <- res_coloc[ , sum(coloc, na.rm = TRUE), by = c("study_exposure", "study_outcome")]
mat_coloc <- format_test_to_matrix(test)
mat_coloc <-mat_coloc[otter_order, otter_order]

test <- res_coloc[ , paste0(sum(coloc, na.rm = TRUE), "/", .N), by = c("study_exposure", "study_outcome")]
mat_proportion <- format_test_to_matrix(test)
mat_proportion <-mat_proportion[otter_order, otter_order]

test <- res_coloc[ , .N, by = c("study_exposure", "study_outcome")]
setnames(test,"N", "V1")
mat_tot <- format_test_to_matrix(test)
mat_tot <-mat_tot[otter_order, otter_order]

coloc_matrix <-  list(mat_coloc = mat_coloc, mat_tot = mat_tot, mat_ratio = mat_ratio, mat_proportion = mat_proportion)

fwrite(res_coloc, "Data/Modified/coloc_res.txt")
fwrite(res_hypr, "Data/Modified/hypr_res.txt")
fwrite(exposure_hypr, "Data/Modified/exposure_hypr.txt")
saveRDS(coloc_matrix, "Data/Modified/coloc_matrix.rds")
message("This script finished without errors")


