#!/usr/bin/env Rscript
library(data.table)
library(GagnonMR)
library(tidyverse)
setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL")

gwasvcf::set_bcftools()
gwasvcf::set_plink()
res_pwas <- fread("Data/Modified/res_pwas.txt")
res_map <- fread("Data/Modified/res_map.txt")
dt_gene_region <- fread("Data/Modified/dt_gene_region.txt")
vep_input <- c(res_pwas$lead_snp.wald, res_map$posprob_coloc.SNP) %>% unique
all_out_vcf <- gwasvcf::query_gwas(vcf = "/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-1-1/trait-1-1.vcf.gz", rsid = vep_input )
                                     
VariantAnnotation::writeVcf(all_out_vcf, file = "Data/Modified/VEPinput.vcf")

system2(command = "sh", args = "Analysis/3c_vep.sh")

###VEPoutput

dt_res <- rbind(res_map, res_pwas, fill = TRUE)
setwd("/home/gagelo01/workspace/Projects/Brain_pQTL")
vepoutput <-fread("Data/Modified/VEPoutput.txt", skip = 39)

{test <- all_out_vcf %>% gwasglue::gwasvcf_to_TwoSampleMR(.) %>% as.data.table
  setnames(vepoutput, "#Uploaded_variation", "Uploaded_variation")
  merge(test[,.(SNP, chr.exposure, pos.exposure)], vepoutput[,.(Uploaded_variation, Location),], 
        by.x = "SNP", by.y = "Uploaded_variation") %>% distinct(.)}
vepoutput <- separate(vepoutput, "Consequence", sep = ",", into = paste0("consequence", 1:4))
vepoutput <- vepoutput[, unique(c(consequence1, consequence2, consequence3, consequence4)), by =c("Uploaded_variation", "Gene")]
setnames(vepoutput, "V1", "consequence")
vepoutput<-vepoutput[!is.na(consequence)]
vepoutput[, consequences := paste(consequence, collapse = ","), by = c("Uploaded_variation", "Gene")]
setnames(vepoutput, "Uploaded_variation", "SNP")
df_VEP <- distinct(vepoutput[, c("SNP", "Gene", "consequences")])
df_VEP <- df_VEP[SNP %in% vep_input, ]
toMatch<- c("missense_variant", "stop_gained", "stop_lost", "start_gained", "start_lost", "frameshift")
df_VEP[ , is_altering_variant := grepl(paste(toMatch,collapse="|"), consequences)]

gene <- c(res_map$hgnc_symbol, res_pwas$hgnc_symbol) %>% unique(.)
ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- biomaRt::getBM(attributes=c("hgnc_symbol","ensembl_gene_id","entrezgene_id"), filters = "hgnc_symbol", mart=ensembl, values=gene)
mirge <- merge(distinct(dt_gene_region[hgnc %in%  gene, .(hgnc, UniProt)]), mapping, by.x = "hgnc", by.y = "hgnc_symbol")

list_res <- map(list(res_pwas, res_map), function(x) {
x <- merge(x, mirge[ensembl_gene_id %in% df_VEP$Gene,.(UniProt, ensembl_gene_id)], by.x = c("UniProt"), by.y = "UniProt",all.x = TRUE)
x <- merge(x, df_VEP, by.x = c("lead_snp.wald", "ensembl_gene_id"), by.y = c("SNP", "Gene"), all.x = TRUE)
x[,ensembl_gene_id := NULL]
x <- distinct(x)
return(x) })
res_pwas <- list_res[[1]]
res_map <- list_res[[2]]

fwrite(res_map, "Data/Modified/res_pwas.txt")
fwrite(res_pwas, "Data/Modified/res_pwas.txt")

message(paste0("This script finished without errors at ", Sys.time()))
