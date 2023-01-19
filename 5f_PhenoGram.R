#!/usr/bin/env Rscript

library(data.table)
library(tidyverse)
  
 setwd( "/mnt/sda/gagelo01/Projects/Brain_pQTL/")
  res_map <- fread("Data/Modified/res_map.txt")
  res_pwas <- fread("Data/Modified/res_pwas.txt")
  dt_res <- fread( "Data/Modified/dt_res_pwas.txt")
  
  # susiesign <- fread( file = "Data/Modified/susiesingstringent.txt" )

 map <- res_map[,.(hgnc_symbol,  study, posprob_colocH4.SNP)]
 map[,phenotype:= "Genome-wide mapping approach"]
 setnames(map, "posprob_colocH4.SNP", "snp")

ntest <- dt_res[!is.na(b.wald), ]$UniProt %>% unique %>% length #number of different proteins tested number of test
pwas <- res_pwas[pval.wald < 0.05/ntest & posprob_coloc_PPH4 > 0.8  & steiger_pval.wald < 0.05 & fstat.exposure.wald > 10, .(hgnc_symbol, study, lead_snp.wald)]
pwas[,phenotype := "Uni-cis PWMR approach"]

# susie<- susiesign[,.(hgnc_symbol, study, susiecoloc.hit_1)]
# setnames(susie, "susiecoloc.hit_1", "snp") 
# susie[,phenotype := "Multi-cis PWMR approach"]

k<-  rbindlist(list(map, pwas), fill = TRUE)
k <- k[, .SD[1], by = c("hgnc_symbol", "phenotype")]      
setnames(k, "hgnc_symbol", "annotation")
k[,study := NULL]
k[is.na(snp), snp := lead_snp.wald]
k[,lead_snp.wald:=NULL]

traduction = fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt")
mirge <- merge(k, traduction[,.(rsid, chr,  position)], by.x = "snp" , by.y = "rsid", all.x = TRUE)  

setnames(mirge, "position", "pos")
mirge[is.na(chr), .N] == 0 #if TRUE excellent
mirge[is.na(pos), .N] == 0 #if TRUE excellent

# mirge2 <- merge(mirge, traduction[,.(rsid, chr,  position)], by.x = c("chr.x", "pos"), by.y = c("chr", "position"), all.x = TRUE)
# mirge2[is.na(snp), snp := rsid]
# mirge2[is.na(rsid), .N == 0] #if TRUE excellent
# setnames(mirge2, "chr.x", "chr")
# mirge2[, c("chr.y", "position", "rsid") := NULL]

fwrite(mirge, "Data/Modified/PhenoGraminput.txt", sep = "\t")


mirge3 <- mirge[, .SD[1,],by = "annotation"]
mirge3[, phenotype := ""]
fwrite(mirge3, "Results/Presentation/PhenoGraminput3.txt", sep = "\t")

# go to http://visualization.ritchielab.org/phenograms/plot
# and add annotation and cytoband

message("This script finished without errors")