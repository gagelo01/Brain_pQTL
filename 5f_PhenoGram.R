#!/usr/bin/env Rscript

library(data.table)
library(tidyverse)
  
 setwd( "/mnt/sda/gagelo01/Projects/Brain_pQTL/")
  res_map <- fread("Data/Modified/res_map.txt")
  res_pwas <- fread("Data/Modified/res_pwas.txt")
  susiesign <- fread( file = "Data/Modified/susiesingstringent.txt" )

 map <- res_map[,.(hgnc_symbol,  study, posprob_colocH4.SNP)]
 map[,phenotype:= "Genome-wide mapping approach"]
 setnames(map, "posprob_colocH4.SNP", "snp")

pwas<- res_pwas[, .(hgnc_symbol, study, chr.exposure, pos.exposure)]
setnames(pwas, c("chr.exposure", "pos.exposure"), c("chr", "pos"))  
pwas[,phenotype := "Uni-cis PWMR approach"]

susie<- susiesign[,.(hgnc_symbol, study, susiecoloc.hit_1)]
setnames(susie, "susiecoloc.hit_1", "snp") 
susie[,phenotype := "Multi-cis PWMR approach"]

k<-  rbindlist(list(map, pwas,susie), fill = TRUE)
k <- k[, .SD[1], by = c("hgnc_symbol", "phenotype")]      
setnames(k, "hgnc_symbol", "annotation")
k[,study := NULL]

traduction = fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt")
mirge <- merge(k, traduction[,.(rsid, chr,  position)], by.x = "snp" , by.y = "rsid", all.x = TRUE)  
mirge[!is.na(snp), c("chr.x","pos") := .(chr.y, position)]    
mirge[is.na(chr.x), .N] == 0 #if TRUE excellent
mirge[is.na(pos), .N] == 0 #if TRUE excellent

mirge2 <- merge(mirge, traduction[,.(rsid, chr,  position)], by.x = c("chr.x", "pos"), by.y = c("chr", "position"), all.x = TRUE)
mirge2[is.na(snp), snp := rsid]
mirge2[is.na(rsid), .N == 0] #if TRUE excellent

setnames(mirge2, "chr.x", "chr")
mirge2[, c("chr.y", "position", "rsid") := NULL]
fwrite(mirge2, "Data/Modified/PhenoGraminput.txt", sep = "\t")


mirge3 <- mirge2[, .SD[1,],by = "annotation"]
mirge3[, phenotype := ""]
fwrite(mirge3, "Results/Presentation/PhenoGraminput3.txt", sep = "\t")

# go to http://visualization.ritchielab.org/phenograms/plot
# and add annotation and cytoband

message("This script finished without errors")