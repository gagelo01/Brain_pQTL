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
dt_res <- rbind(res_map, res_pwas, fill = TRUE) %>% distinct
vepoutput <-fread("Data/Modified/VEPoutput.txt", skip = 39)
setnames(vepoutput, c("#Uploaded_variation", "Consequence"), c("Uploaded_variation", "consequence"))
clean_vepoutput<-function(vepoutput){
  k<-vepoutput[,.(Uploaded_variation, Gene, consequence)]
  nom_var<-"consequence"
  num<-k[,str_count(get(nom_var), pattern = ",")%>%max+1]
  k <- separate(data = k, col = nom_var, sep = ",", into = paste0(nom_var, 1:num))
  k <- melt(k, id.vars = c("Uploaded_variation","Gene"), measure.vars = paste0(nom_var, 1:num))
  k <- distinct(k[!is.na(value),.(Uploaded_variation, Gene, value)])
  my_aggregate_fun<-function(x){paste(x, collapse = ", ")}
  k[,consequence:="consequence"]
  k <- dcast(k, Uploaded_variation + Gene ~ consequence, value.var = "value", fun.aggregate = my_aggregate_fun)
  toMatch<- c("missense_variant", "stop_gained", "stop_lost", "start_gained", "start_lost", "frameshift")
  k[ , is_altering_variant := grepl(paste(toMatch,collapse="|"), consequence)]
  ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mapping <- biomaRt::getBM(attributes=c("hgnc_symbol","ensembl_gene_id"), filters = "ensembl_gene_id",
                            mart=ensembl, values= unique(k$Gene))
  k <- merge(k, mapping, by.x = "Gene",  by.y = "ensembl_gene_id", all.x = TRUE)
  setnames(k ,"hgnc_symbol", "hgnc")
  return(k)
}

vepoutput_clean <- clean_vepoutput(vepoutput = vepoutput)

fwrite(vepoutput_clean, "Data/Modified/vepoutput_clean.txt")

message(paste0("This script finished without errors at ", Sys.time()))
