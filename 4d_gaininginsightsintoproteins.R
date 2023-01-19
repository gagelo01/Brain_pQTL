#!/usr/bin/env Rscript
library(tidyverse)
library(data.table)
library(GagnonMR)

setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL")
res_specificity <- fread("Data/Modified/res_specificity.txt")
res_pwas <- fread( "Data/Modified/res_pwas.txt")
dt_res <- fread("Data/Modified/dt_res_pwas.txt")
ntest <- dt_res[!is.na(b.wald), ]$UniProt %>% unique %>% length #number of different proteins tested number of test
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")

####
k <- rbindlist(list(res_pwas, res_specificity), fill = TRUE)
kcause <- k[pval.wald < 0.05/ntest & posprob_coloc_PPH4 > 0.8  & steiger_pval.wald < 0.05 & fstat.exposure.wald > 10, unique(UniProt), by = "study"]
listInput <- lapply(split(kcause, kcause$study), function(x) x$V1)
#########
myfromList <- function (input) {
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  data$ELEMENTS <- elements
  return(as.data.table(data))
}

test<- myfromList(listInput)
colini<-names(listInput)
test[,(colini):=lapply(.SD, as.numeric), .SDcols = colini]
test[banner == 1, ELEMENTS[apply(.SD, 1, sum)==1], .SDcols = colini]
k <- test[banner == 1 & rosmap == 1, ELEMENTS[apply(.SD, 1, sum)==2], .SDcols = colini]
dt_gene_region[UniProt %in% k, unique(hgnc)]


#Genes known to be involved in obesity¸
geneset1 <- c("LEP", "LEPR", "POMC", "MC4R", "SIM1", "NTRK2", "KSR2", "CPE", "PCSK1", "BDNF", "SH2B1", "TUB") #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6226269/
geneset2 <- c("SLC6A17", "RAPGEF3", "PRKAG1", "RAB21", "KSR2", "MAP1A", "MC4R", "GIPR", "AKDH3A1", "ZFR2", "ACHE", "ANGPTL7", "ZNF169") #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5945951/#SD2
geneset3 <- c("LEP", "LEPR", "POMC", "BDNF", "MC4R", "PCSK1", "SIM1", "BDNF", "FTO", "NRXN3", "NPC1", "NEGR1", "MTCH2", "GNPDA2", "APOE", "CD38", "SIRT1", "TNFA", "SERPINB2", "TREM2", "SYT4", "FMR1", "TET3",
              "MC4R", "POMC", "FTO", "NRXN3", "NPC1", "NEGR1", "GNPDA2", "MTCH2", "ETV5", "NEGR1", "NRXN3", "CADM2", "GRID1", "ELAVL4", "SCG3", "ETV5", "HNF4G", "TLR4", "ADCY3") #https://www.frontiersin.org/articles/10.3389/fnins.2020.00863/full
geneset4<-c("ADCY3", "MYT1L", "POU3F2", "GRPR", "LRP2", "ADCY3", "BDNF", "CPE", "GRPR", "LEP", "LEPR", "LRP2", "MC3R", "MC4R", "MRAP2", "MYT1L", "NPY", "NTRK2", "PCSK1", "POMC", "SH2B1", "SIM1", "TUB", "CPE", "GRPR", "LEP", "LRP2", "MC3R", 
            "MRAP2", "MYT1L", "NPY", "NTRK2", "PCSK1", "POMC", "SH2B1", "SIM1", "TUB") #https://www.frontiersin.org/articles/10.3389/fendo.2020.00081/full
genetable <- fread("Data/Raw/geneobesitytablemid28910990.txt")
geneset <- c(geneset4, geneset2, geneset3, geneset1, genetable$`Gene symbol`) %>% unique

saveRDS(geneset, file = "Data/Modified/known_geneset.rds")
protset <- dt_gene_region[hgnc%in%geneset, unique(UniProt)]
protset <- test[ELEMENTS %in% protset, ]$ELEMENTS
dt_gene_region[UniProt%in%protset,unique(hgnc)]
