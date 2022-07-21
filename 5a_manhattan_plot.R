#!/usr/bin/env Rscript
library(CMplot)
library(data.table)
library(tidyverse)
setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL")
source("Analysis/CMplot_eloi.R")
res_map <- fread( "Data/Modified/res_map.txt")
tomanhatthan<- res_map[ , .SD[which.max(posprob_coloc_PPH4)] , by = "hgnc_symbol", .SDcols = c("posprob_colocH4.SNP")]

outcome_bmi_ukb_giant_full <- fread( "Data/Modified/outcome_bmi_ukb_giant_full.txt")
data <- outcome_bmi_ukb_giant_full[!(is.na(chr.outcome) | is.na(pos.outcome) | is.na(pval.outcome)) & pval.outcome < 5e-5,]
data <- data[,.(SNP, chr.outcome, pos.outcome, pval.outcome)]

setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL/Results")
CMplot_eloi(as.data.frame(data), type="p", plot.type="m", LOG10=TRUE, threshold=5e-08, amplify=FALSE, memo="", dpi=300, verbose=TRUE, width=14,
            height=6, col=c("grey70", "grey90"), threshold.lwd=2, cex=0.4, highlight.cex =  0.6, highlight=tomanhatthan$posprob_colocH4.SNP,
            highlight.text=tomanhatthan$hgnc_symbol, highlight.col = "red", 
            highlight.text.xadj =rep(-1, length(tomanhatthan$hgnc_symbol)), highlight.text.yadj =rep(1, length(tomanhatthan$hgnc_symbol)),
            file.output = TRUE, arrow_col_eloi = "black", highlight.text.cex=0.7)

setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL/Results/Presentation")
CMplot_eloi(as.data.frame(data), type="p", plot.type="m", LOG10=TRUE, threshold=5e-08, amplify=FALSE, memo="", dpi=300, verbose=TRUE, width=14,
            height=6, col=c("grey70", "grey90") , threshold.lwd=2, cex=0.4,file.output = TRUE, arrow_col_eloi = "black")

message("This script finished without errors")

