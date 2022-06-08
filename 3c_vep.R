library(data.table)
library(GagnonMR)
library(tidyverse)
setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL")

gwasvcf::set_bcftools()
gwasvcf::set_plink()
res_pwas <- fread("Data/Modified/res_pwas.txt")
res_map <- fread("Data/Modified/res_map.txt")
vep_input <- c(res_pwas$lead_snp.wald, res_map$posprob_coloc.SNP) %>% unique
all_out_vcf <- gwasvcf::query_gwas(vcf = "/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-1-1/trait-1-1.vcf.gz", rsid = vep_input[1] )
                                     
VariantAnnotation::writeVcf(all_out_vcf, file = "Data/Modified/VEPinput.vcf")
fwrite(list(vep_input), "Data/Modified/VEPinput.vcf")

system2(command = "sh", args = "Analysis/3c_vep.sh")

# then go see this website https://www.ensembl.org/info/docs/tools/vep/index.html
#follow everything with default settings except change for a 200Kb window


###VEPoutput
###VEPoutput

dt_res <- rbind(res_map, res_pwas, fill = TRUE)
setwd("/home/gagelo01/workspace/Projects/Brain_pQTL")
vepoutput <-fread("Data/Modified/VEPoutput.txt", skip = 39)

vepoutput<- vepoutput[SYMBOL %in% dt_res$exposure, ]

vepoutput <- separate(vepoutput, "Consequence", sep = ",", into = paste0("consequence", 1:4))

vepoutput <- vepoutput[, unique(c(consequence1, consequence2, consequence3, consequence4)), by =c("#Uploaded_variation", "SYMBOL")]

setnames(vepoutput, "V1", "consequence")
vepoutput<-vepoutput[!is.na(consequence)]

vepoutput[, consequences := paste(consequence, collapse = ","), by = c("#Uploaded_variation", "SYMBOL")]
setnames(vepoutput, "#Uploaded_variation", "SNP")

vepoutput[,SNP_exposure := paste0(SNP, "-", SYMBOL  )]
grepl("-",dt_res$exposure) %>% any # if FALSE means I do not need to grep, buy can do exact matching

res_pwas[, SNP_exposure := paste0(lead_snp.wald, "-", exposure)]
res_map[, SNP_exposure := paste0( posprob_coloc.SNP, "-", exposure)]

df_VEP <- distinct(vepoutput[, c("SNP", "SYMBOL", "SNP_exposure","consequences")])
toMatch<- c("missense_variant", "stop_gained", "stop_lost", "start_gained", "start_lost", "frameshift")
df_VEP[ , is_altering_variant := grepl(paste(toMatch,collapse="|"), consequences)]

df_VEP[, .(SYMBOL, is_altering_variant)][is_altering_variant == TRUE]
res_pwas <- merge(res_pwas, df_VEP[,!c("SNP", "SYMBOL")], by = "SNP_exposure")

res_map <- merge(res_map, df_VEP[,!c("SNP", "SYMBOL")], by = "SNP_exposure")
res_pwas[is_altering_variant == TRUE]
res_map[is_altering_variant == TRUE]

fwrite(res_map, "Data/Modified/res_pwas.txt")
fwrite(res_pwas, "Data/Modified/res_pwas.txt")

message(paste0("This script finished without errors at ", Sys.time()))