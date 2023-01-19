#!/usr/bin/env Rscript
suppressWarnings(suppressPackageStartupMessages({
library(data.table)
  library(tidyverse)
  library(GagnonMR)
  library("xlsx")

}))
setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
res_map <- fread("Data/Modified/res_map.txt")
res_pwas <- fread("Data/Modified/res_pwas.txt")
BMI_loci <- fread("Data/Modified/BMI_loci.txt")
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")
df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
dt_res <- fread( "Data/Modified/dt_res_pwas.txt")
posprob_matrix <- readRDS("Data/Modified/coloc_matrix.rds")
hypr_res <- fread( "Data/Modified/hypr_res.txt")
hypr_res<- distinct(hypr_res)
res_full <- rbindlist(list(res_pwas, res_map), fill = TRUE) %>% distinct(.)
dt_replication_brainpqtl <- fread( "Data/Modified/dt_replication_brainpqtl.txt")
res_specificity <- fread("Data/Modified/res_specificity.txt")
coloc_matrix <- readRDS( "Data/Modified/coloc_matrix.rds")
known_geneset <- readRDS("Data/Modified/known_geneset.rds")
vepoutput <- fread("Data/Modified/vepoutput_clean.txt")
data_united_tpm <- fread("Data/Modified/data_united_tpm.txt")
df_index <- data.table::fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")

# res_yang <- fread( "Data/Modified/res_yang.txt")
# dt_res_hypr <- fread( "Data/Modified/dt_hypr_dietaryhabits.txt")
wingobmipwas <- xlsx::read.xlsx( "/mnt/sda/gagelo01/Projects/Brain_pQTL/Data/Raw/41593_2021_832_MOESM2_ESM.xlsx", 
                                 sheetName ="SupTb11", startRow = 4) %>% as.data.table
compare <- fread( "Data/Modified/compare_tpm.txt")
data_united <- fread( "Data/Modified/data_united_tpm.txt")
# grs <- fread( "Data/Modified/grs_results_norm.txt")

return_format_data<-function(data) {
  return(data[, paste0(round(b, digits = 2), " 95% CI=", round(lci, digits = 2), "-",  round(uci, digits = 2), ", p=",pval %>% formatC(., format = "e", digits = 1))])
}


#NUmber of genes in the vicinity of BMI loci
gene_region <- dt_gene_region[, .(hgnc, gene_region)] %>% distinct
gene_region <- separate(gene_region, col = "gene_region", into = c("chr", "start", "end"))
gene_region[, c("chr", "start", "end") := lapply(.SD, as.numeric), .SDcols = c("chr", "start", "end")]

BMI_loci_split <- split(BMI_loci, by = "SNP")
Ngene <- lapply(BMI_loci_split, function(bon) {
  Ngene <- gene_region[chr == bon$chr.outcome & start < bon$pos.outcome & bon$pos.outcome < end, length(unique(hgnc))]
  return(data.frame(chr.outcome = bon$chr.outcome, pos.outcome = bon$pos.outcome, Ngene = Ngene))
}) %>% rbindlist
BMI_loci <- merge(BMI_loci, Ngene, by = c("chr.outcome", "pos.outcome"))
tablelocigene <- BMI_loci$Ngene %>% table ###nloci

#Figure 1
dt_res[, length(unique(hgnc_symbol)), by = "study"]
dt_res$hgnc_symbol %>% unique(.) %>% length
#Abstract
BMI_loci[,.N]
res_map$hgnc_symbol %>% unique %>% length
setdiff(c(c(res_pwas$UniProt) %>% unique), res_map$UniProt) %>% length
paste0(hypr_res[traits!="None"&grepl("GTEX",traits),gene_investigated%>%unique%>% length ], "/",
hypr_res[,gene_investigated%>%unique%>% length ])

c(res_map$UniProt, res_pwas$UniProt) %>% unique %>% length

###INtroduction
dt_gene_region$UniProt %>% unique %>% length

###Results
####MAP
####para1
BMI_loci[,.N]
df_index[id=="trait-1-1", sample_size]
dt_gene_region$UniProt %>% unique %>% length
res_map$hgnc_symbol %>% unique %>% length

###para 2
{res_map[ posprob_colocH4.SNPexplained_var > 0.6 ,] 
tomanhatthan<- res_map[ , .SD[which.max(posprob_coloc_PPH4)] , by = "hgnc_symbol", .SDcols = c("posprob_colocH4.SNP")]
ldmat <- ieugwasr::ld_matrix_local(tomanhatthan$posprob_colocH4.SNP, plink_bin = genetics.binaRies::get_plink_binary(), 
                                   bfile = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs")
ldmat <- ldmat^2 
diag(ldmat)<-0
max(ldmat) #
nrow(ldmat)} # It never returned more than one gene

dt_res[!is.na(b.wald), length(unique(UniProt))]
ntest <- dt_res[!is.na(b.wald), ]$UniProt %>% unique %>% length #number of different proteins tested number of test
res_pwas[,length(unique(UniProt))]
0.05/ntest
res_pwas[,any(steiger_pval.wald),by = "UniProt"][,.N]
res_pwas[fstat.exposure.wald > 10, .N]
res_pwas[is_in_pleiotropic_region == TRUE,]
vep_output_clean <- merge(vepoutput[,.(Uploaded_variation, hgnc, consequence, is_altering_variant)], 
                          res_pwas[,.(lead_snp.wald, hgnc_symbol, study)],
                          by.x = c("Uploaded_variation", "hgnc"),
                          by.y = c("lead_snp.wald", "hgnc_symbol"))
vep_output_clean[is_altering_variant == TRUE,]$hgnc %>% unique %>% length

dt_replication_brainpqtl$category %>% table

#Para 3
paste0(data_united_tpm$Description %>% unique(.) %>% length(.), "/",
      length(unique(c(res_pwas$hgnc_symbol, res_map$hgnc_symbol))))
data_united_tpm[fTau<0.7, length(unique(Description))]
{data_united_tpm[, in_brain := grepl("brain", tolower(tissue))]
k<-data_united_tpm[, mean(Averaged.TPM), by = c("Description", "in_brain")]
k<-dcast(k, Description ~ in_brain, value.var = "V1")
setnames(k, c("TRUE","FALSE"), c("in_brain", "not_in_brain"))
k[Description%in%data_united_tpm[fTau>0.7,unique(Description)] & in_brain > not_in_brain,]$Description}
#para4
res_specificity[UniProt%in%res_pwas$UniProt & study %in% c("deCODE","FENLAND"), length(unique(UniProt))]
res_specificity[UniProt%in%res_pwas$UniProt & study %in% c("deCODE","FENLAND")&posprob_coloc_PPH4>0.8&pval.wald<0.05/ntest, length(unique(UniProt))]
res_specificity[UniProt%in%res_pwas$UniPro & study%in% "eQTLGen", length(unique(UniProt))]
res_specificity[UniProt%in%res_pwas$UniPro & study%in% "eQTLGen"&posprob_coloc_PPH4>0.8&pval.wald<0.05/ntest, length(unique(UniProt))]

#Para 5
dt_res[!is.na(b.wald),length(unique(UniProt))]
coloc_matrix$mat_ratio
df_index[consortium == "eQTLGen", length(unique(trait))]
df_index[consortium == "FENLAND", .N]
df_index[consortium == "deCODE", length(unique(note))]
df_index[consortium == "deCODE", ]

#Para 5 
hypr_res[traits!="None"&grepl("GTEX",traits),gene_investigated%>%unique%>% length ]
hypr_res[,gene_investigated%>%unique%>% length ]
hypr_res[gene_investigated == "ADCY3", ]$regional_prob

#discussion
#para 1
c(res_pwas$UniProt, res_map$UniProt) %>% unique %>% length
setdiff(unique(susiesign$exposure), res_pwas$exposure) %>% length

# para 2
c("DOC2A", "COMT", "ADCY3", "CAMKK2") %in% c(res_pwas$hgnc_symbol, res_map$hgnc_symbol)
intersect(known_geneset, c(res_pwas$hgnc_symbol, res_map$hgnc_symbol))

# para 3
coloc_matrix$mat_ratio

#para 4 comparison with wingo
wingobmipwas <- wingobmipwas[!is.na(Gene),]
wingobmipwas$Gene %>% unique(.) %>% length
listprotein <- unique(c(res_pwas$hgnc_symbol, res_map$hgnc_symbol))                
intersect(listprotein, wingobmipwas$Gene) %>% unique %>% length
setdiff(listprotein, wingobmipwas$Gene  )  %>% unique %>% length 
dt_res[hgnc_symbol %in% unique(wingobmipwas$Gene),][posprob_coloc_PPH4 >0.8, length(unique(hgnc_symbol))]

#para 5 limits
dt_gene_region$UniProt %>% unique %>% length
dt_res[pval.exposure.wald < 5e-8, length(unique(UniProt))]

# para 6
c(res_map$UniProt, res_pwas$UniProt)  %>% unique %>% length
