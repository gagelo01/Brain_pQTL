#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(GagnonMR)
setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL")

tissue <- fread("Data/Modified/MAGMA_tissue_res.txt", skip = 5)
# tissue_toremove<- c("Vagina", "Uterus", "Testis", "Prostate", "Ovary", "Breast_Mammary_Tissue", "Fallopian_Tube", "Cervix_Ectocervix",
#                     "Cervix_Endocervix", "Bladder", "Kidney_Cortex", "Kidney_Medulla")
# 
# all(tissue_toremove %in% M$FULL_NAME) #if TRUE good
# M2 <- tissue[!(FULL_NAME %in% tissue_toremove), ]

celltype <- fread("Data/Modified/magma_celltype_GSE67835_Human_Cortex.gsa.out", skip = 5)
setnames(celltype, "VARIABLE", "FULL_NAME")

histogramm_tissue_celltype <- function(M) {
M[,logp := -log10(P)]
M[, tissue_name := gsub("_", " ", FULL_NAME)]
M <- M[order(-logp)]
M[,tissue_name := factor(tissue_name, levels = unique(tissue_name))]
bonfp <- M[,-log10(0.05/.N)] #bonferroni correction threshold
nomp <- -log10(0.05)
M[, colourbar := logp %>% ifelse(. > bonfp, "red3", .) %>% ifelse(logp > nomp & logp < bonfp, "blue", .) %>% ifelse(. != "red3" & . != "blue", "grey", .)]
M[,colourbar := factor(colourbar, levels = c("red3", "blue", "grey"))]

ggplot(M, aes(x = tissue_name, y =logp, fill = colourbar)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(face = "bold", angle = 85, hjust = 1)) +
  scale_fill_manual("legend", values = levels(M$colourbar)[levels(M$colourbar) %in% unique(M$colourbar)]) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  geom_hline(yintercept = c(bonfp, nomp), linetype = "dashed", color = c("red3", "blue")) + 
  guides(fill="none") +
  ylab("-log 10 P-Value") +
  xlab("")

}

histogramm_tissue_celltype(M = tissue)

ggsave("Results/barplot_tissular_specificiy_.png", 
       width = 831/72,
       height = 466/72,
       units = "in")

####Cell type

histogramm_tissue_celltype(M = celltype)

ggsave("Results/barplot_celltype_enrichment.png",
                           width=290/72,
                           height=250/72,
                           units = "in")


