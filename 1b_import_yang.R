#!/usr/bin/env Rscript
#import instrument
library(data.table)
library(readxl)
library(TwoSampleMR)
library(tidyverse)
library(AbHAC)

setwd("/mnt/sde/gagelo01/Projects/Brain_pQTL")

yang<-read_xlsx("Data/Raw/f3MultitissuepQTLsupTablesforNatureNeuroAug312020.xlsx",
                sheet = 6 , range = "A2:O2680")
setDT(yang)
#rsids

snps <- yang[,"start"]
colnames(snps) <- "pos"	# Car on veut que zgrep trouve tous les items de la liste DONT l'entête qui contient colname_snps...
fwrite(snps, file = "snps_list_part1.txt")		# writing de la liste en format .txt in the current directory

all_out <-  fread(cmd = paste0("zgrep -w -f snps_list_part1.txt ", "/mnt/sde/couchr02/rsids/rsid_pos_86millions.txt"))
all_out[, chr_pos := paste0(chr, ":", pos)]

all(yang$start == yang$end)
yang[, chr_pos := paste0(chr, ":", start)]
all_out <- all_out[chr_pos %in% yang$chr_pos]

mirge <- merge(yang, all_out[,c(4,1) ], by.x = "chr_pos", by.y = "chr_pos", all.x = TRUE)

mirge<-mirge[!is.na(rsid),]

if (file.exists("snps_list_part1.txt")) {
  file.remove("snps_list_part1.txt")
}

##get uniprot id
dt_aptamers <- as.data.table(readat::aptamers)
mirge <- merge(mirge, distinct(dt_aptamers[, c("Target", "UniProt")]), by.x = "TargetName", by.y = "Target", all.y = FALSE, all.x = FALSE)

mirge[, pQTL_affectGene := AbHAC::uniprot.to.hgnc(UniProt)]
mirge[,units := "SD"]
mirge[,z := abs(qnorm(p_value))] #z-score from one tailed p-value
mirge[,se := abs(BETA/z)]
mirge[, samplesize := 459] #We generated data for 1,305 proteins using an aptamer-based approac brain (n=459) samples 
mirge[,id := UniProt]

#add eaf #genomic positions (start and end coordinates) was performed using gencode v30 liftover to hg19/GRCh37.
gencode <- data.table::fread("/home/couchr02/Mendel_Commun/Nicolas/GTEx/gencode.v19.genes.v7.patched_contigs.txt")
traduction <-  fread("/mnt/sde/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt")
setorder(traduction, chr, position)
traduction <- traduction[, .(rsid, chr, position, a0, a1, EUR)]
traduction[ , EUR := EUR %>% ifelse(. == 0, 0.001, .) %>% ifelse(. == 1, 0.999, .)]

window <- 2e5
gene_coords <- gencode
gene_coords[, start := min(start), by = "gene_name"]
gene_coords[, end := max(end), by = "gene_name"]
gene_coords <- distinct(gene_coords)
gene_coords[, start_cis := (start - window/2) %>% ifelse(.<1, 1, .) ]
gene_coords[, end_cis := end + window/2]

traduction[,chr:=as.integer(chr)]
traduction[,position:=as.integer(position)]
mirge2 <- merge(mirge, traduction[,-c("rsid")], by.x = c("chr", "start"), by.y = c("chr", "position"))
mirge2 <- mirge2[(altAllele == a0 | altAllele == a1) & (refAllele == a0 | refAllele == a1) & a0 != a1 & refAllele != altAllele,]
mirge2 <- mirge2[chr %in% 1:22,]
mirge2[altAllele == a0, BETA := BETA*-1]
mirge2[altAllele == a0, altAllele := a1]
mirge2[refAllele == a1, refAllele := a0]
mirge3 <- merge(mirge2, gene_coords[,-c("chr", "start", "end")], by.x = "pQTL_affectGene", by.y = "gene_name")
mirge3[, pqtl_type := ifelse( (start > start_cis) & (start < end_cis), "cis", "trans")]
mirge3[, eaf := ifelse(EUR < 0.5, MAF, 1 - MAF)]
#format_data
mirge3 <- distinct(mirge3)
inst_yang <- format_data(mirge3,
                         type = "exposure",
                         phenotype_col = "UniProt" ,
                         snp_col = "rsid",
                         beta_col = "BETA",
                         se_col = "se",
                         effect_allele_col = "altAllele",
                         other_allele_col = "refAllele",
                         samplesize_col = "samplesize",
                         pval_col = "p_value",
                         units_col = "units",
                         id_col = "id",
                         gene_col = "pQTL_affectGene",
                         chr_col =  "chr",
                         pos_col = "start",
                         eaf_col = "eaf")

inst_yang <- merge(inst_yang, mirge3[, c("rsid", "UniProt", "pqtl_type", "MAF" )], by.x = c("SNP", "exposure"), by.y = c("rsid", "UniProt"))
setDT(inst_yang)

inst_yang[, exposure := paste0("yang", "_", exposure)]
inst_yang<-separate(inst_yang, exposure, into = c("study", "protein"), sep = "_", remove = FALSE)

inst_yang<- GagnonMR::remove_gene_region(inst_yang)
inst_yang[, MAF := ifelse(eaf.exposure<0.5, eaf.exposure, 1 - eaf.exposure)]
inst_yang[,max(pval.exposure)]
inst_yang[, id.exposure := exposure]

inst_clump <- GagnonMR::clump_data_local(inst_yang)

fwrite(inst_clump, "Data/Modified/inst_yang_clump.txt")

message("This script finished withour errors")


