#!/usr/bin/env Rscript
#import instrument
library(data.table)
library(readxl)
library(TwoSampleMR)
library(biomaRt)
library(GagnonMR)
library(furrr)
library(tictoc)
setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL")
rosmap <- fread("Data/Raw/Robins_DLPFC_pQTLs/ROSMAP_DLPFC_pQTLs.csv")
banner <- fread("Data/Raw/Robins_DLPFC_pQTLs/BannerBBDP_DLPFC_pQTLs.csv")
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

#uniprot in gene
uniprot_dictionnary <- fread("/mnt/sda/boujer01/Drug_Targets/Uniprot_ID/uniprot_final.txt")
colnames(uniprot_dictionnary)<- c("uniprot", "gene_name", "n_aa")

uniprot_vec <- unique(c(banner$UNIPROT, rosmap$UNIPROT))
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

uniprot_filters <- c("uniprotswissprot", "uniprot_gn_id", "uniprot_gn_symbol", "uniprot_isoform",  "uniprotsptrembl")
uniprot_biomart <- data.table(hgnc_symbol = "hgnc_symbol", uniprot_gn_id = "uniprot_gn_id")
for(i in 1:length(uniprot_filters)) {
  values <- uniprot_vec[!(uniprot_vec %in% uniprot_biomart$uniprot_gn_id)]
  k <- getBM(attributes=c(uniprot_filters[i], "hgnc_symbol"),
             filters =  uniprot_filters[i] , values=values, mart=ensembl)
  setDT(k)
  setnames(k, uniprot_filters[i], "uniprot_gn_id")
  uniprot_biomart <- rbind(uniprot_biomart, k)
}
values <- uniprot_vec[!(uniprot_vec %in% uniprot_biomart$uniprot_gn_id)]

vec_res <- vector(mode = "list", length = length(values))
for(i in 1:length(values)) {
  x  <- uniprot_dictionnary[grepl(values[i],uniprot),] #gsub("-.*", "",values[i])
  x[,uniprot_gn_id := values[i]]
  setnames(x, "gene_name",c("hgnc_symbol"))
  vec_res[[i]] <- x[, .(hgnc_symbol, uniprot_gn_id)]
}
dt_res <- rbindlist(vec_res)
uniprot_biomart <- rbind(uniprot_biomart, dt_res)
uniprot_biomart <- uniprot_biomart[!(hgnc_symbol == "hgnc_symbol"),]
uniprot_biomart <- uniprot_biomart[hgnc_symbol != "",]
uniprot_vec[!(uniprot_vec %in% uniprot_biomart$uniprot_gn_id)]
double_hgnc <-uniprot_biomart[, .N, by = "uniprot_gn_id"][N>1]$uniprot_gn_id
ubd <- uniprot_biomart[uniprot_gn_id %in% double_hgnc, ][, .SD[2],by = "uniprot_gn_id"]
uniprot_biomart <- uniprot_biomart[!(uniprot_gn_id %in% ubd$uniprot_gn_id),] 
uniprot_biomart <- rbind( uniprot_biomart, ubd)
rosmap[!(UNIPROT %in% uniprot_biomart$uniprot_gn_id), ]
banner[!(UNIPROT %in% uniprot_biomart$uniprot_gn_id), unique(UNIPROT)]

all(rosmap$UNIPROT %in% uniprot_biomart$uniprot_gn_id )
all(uniprot_biomart$uniprot_gn_id %in% banner$UNIPROT)
# obtain_genecard_from_uniprot <- function(uniprot_vec, mart_human) {
# uniprot_vec <- gsub("-.", uniprot_vec)
# gene_coords_mart <- getBM(attributes=c("hgnc_symbol", "uniprot_gn_id"),
#                           filters="uniprot_gn_id", values=uniprot_vec, mart=mart_human)
# 
# setDT(gene_coords_mart)
# gene_coords_mart <- gene_coords_mart[hgnc_symbol != "",]
# uniprotendouble<- gene_coords_mart[, length(unique(hgnc_symbol)), by = "uniprot_gn_id"][V1>1,]$uniprot_gn_id
# if(length(uniprotendouble)>0) {
# gene_coords_mart[uniprot_gn_id %in% uniprotendouble, ]
# index <- sapply(as.list(gene_coords_mart[uniprot_gn_id %in% uniprotendouble, ]$hgnc_symbol), function(x) 
#   !(from_genecard_to_generegion(genecard_name = x) %>% is.null))
# gene_coords_mart <- gene_coords_mart[!(hgnc_symbol %in% c("BCL2L2", "AARSD1", gene_coords_mart[uniprot_gn_id %in% uniprotendouble, ][!index,]$hgnc_symbol)),]
# }
# gene_coords_mart[,hgnc_symbol := paste(unique(hgnc_symbol), collapse = "-"),by = "uniprot_gn_id"]
# gene_coords_mart <- distinct(gene_coords_mart)
# double_uniprot <- gene_coords_mart[, .N, by = "uniprot_gn_id"][N>1]
# stopifnot(nrow(gene_coords_mart[uniprot_gn_id %in% double_uniprot$uniprot_gn_id, ][order(uniprot_gn_id)]) == 0)#SHOULD BE TRUE
# return(gene_coords_mart)
# }

# k_rosmap <- obtain_genecard_from_uniprot(uniprot_vec = unique(rosmap$UNIPROT), mart_human = human)
# k_banner <- obtain_genecard_from_uniprot(uniprot_vec = unique(banner$UNIPROT), mart_human = human)
# vec_uniprot <- rosmap[!(UNIPROT %in% k_rosmap$uniprot_gn_id),unique(UNIPROT)]
# vec_res <- vector(mode = "list", length = length(vec_uniprot))
# for(i in 1:length(vec_uniprot)) {
#   vec_res[i]  <- uniprot_dictionnary[grepl(vec_uniprot[i],uniprot),]
#   
# }
# 
# fwrite(k_rosmap, "Data/Modified/k_rosmap.txt")
# fwrite(k_banner, "Data/Modified/k_banner.txt")


#align the whole dataset with traduction
traduction <-  fread("/mnt/sda/couchr02/1000G_Phase3/1000G_Phase3_b37_rsid_maf.txt")
setorder(traduction, chr, position)
traduction[ , EUR := EUR %>% ifelse(. == 0, 0.001, .) %>% ifelse(. == 1, 0.999, .)]

harmonisewithtrad <- function(all_out, traduction) {
all_out <- merge(all_out, traduction, by.x = c("CHR", "POS"), by.y = c("chr", "position"), all = FALSE) 
all_out <- all_out[(ALT == a0 | ALT == a1) & (REF == a0 | REF == a1) & a0 != a1  & ALT != REF, ] #because low number removed, coded on the forward strand
all_out <- all_out[CHR %in% 1:22, ]
all_out[, CHR := as.integer(CHR)]
all_out[ALT == a0, BETA := BETA*-1]
all_out[ALT == a0, ALT := a1]
all_out[REF == a1, REF := a0] #less than mrbase, possibly because I do not use the same traduction file
return(all_out)
}

rosmap<- harmonisewithtrad(all_out = rosmap, traduction = traduction)
rosmap[,study := "rosmap"]
split_rosmap <- split(rosmap, rosmap$UNIPROT)

banner<- harmonisewithtrad(all_out = banner, traduction = traduction)
banner[,study := "banner"]
split_banner <- split(banner, banner$UNIPROT)

# create_newk <- function(pqtl, k) {
  # nsnp<- pqtl[, .N, by = "UNIPROT"]
#   k <- merge(k, nsnp, by.x = "uniprot_gn_id", by.y = "UNIPROT")
#   nsample<-pqtl[, max(unique(N)), by = "UNIPROT"]
#   k <- merge(k, nsample, by.x = "uniprot_gn_id", by.y = "UNIPROT")
#   k[,study := pqtl$study[1]]
#   pqtl_copy  <- pqtl[,unique(paste0(CHR, ":", min(POS), "-", max(POS))), by = "UNIPROT"]
#   setnames(pqtl_copy, "V1", "chrompos")
#   k <- merge(k, pqtl_copy, by.x = "uniprot_gn_id", by.y = "UNIPROT")
#   return(k)
# }

# k_rosmap <- create_newk(pqtl = rosmap, k = k_rosmap)
# k_banner <- create_newk(pqtl = banner, k = k_banner)


create_newrow <- function(k, id) {
 k[, chrpos := paste0(unique(CHR), ":", min(POS), "-", max(POS)) , by = "UNIPROT"]
 k[,nsnp:=  .N, by = "UNIPROT"]
 k[, nsample := max(N), by = "UNIPROT"]
 k<- k[,.SD[1],by ="UNIPROT"]
 k <- merge(k, uniprot_biomart, by.x = "UNIPROT", by.y = "uniprot_gn_id", all.x =  TRUE)
 k[is.na(hgnc_symbol), hgnc_symbol := "nonavailable"]
newrow<- data.frame( id  = paste0(id, 1:(k[, length(unique(UNIPROT))] )), trait = k$UNIPROT, group_name = "None", year = 2021,
                       author = "Robins Chloe",
                       consortium =    k[1,study], sex = "Males and Females", population = "European", unit = "SD",
                       nsnp = k$N,sample_size = k$nsample, initial_build = "HG19/GRCh37", category = "protein", pmid = 33571421,
                       sd = 1, note = k[,paste0(hgnc_symbol, "_", chrpos )], ncase = NA, ncontrol = NA)
  return(newrow)
}

df_index <- fread("/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")
df_index <- df_index[pmid != 33571421,]

newrow <- create_newrow( k = rosmap, id = "prot-2-")
df_index <- rbind(df_index, newrow)
newrow <- create_newrow( k = banner, id = "prot-3-")
df_index <- rbind(df_index, newrow)
setDT(df_index)

pqtl <- split_rosmap[[2]]
brain_pqtl_formatvcf_wrapper <- function(pqtl, df_index) {
  
  message(paste0("Initialising ", pqtl$UNIPROT[1]))
  
  formattovcf_createindex(
    all_out = pqtl,
    snp_col = "rsid",
    outcome_name = pqtl[1,]$UNIPROT,
    beta_col = "BETA",
    se_col = "SE",
    pval_col =  "P",
    eaf_col = "EUR",
    effect_allele_col = "ALT",
    other_allele_col = "REF",
    ncase_col =  NULL,
    ncontrol_col =  NULL,
    samplesize_col = "N",
    chr_col = "CHR",
    pos_col = "POS",
    units = "SD",
    traduction = NULL,
    df_index = df_index,
    group_name = "None",
    year = 2021,
    author = "Robins Chloe",
    consortium = pqtl$study[1],
    sex = "Males and Females",
    population = "European", 
    initial_build = "HG19/GRCh37",
    category = "protein",
    pmid = 33571421, 
    note = df_index[consortium == pqtl$study[1] & trait == pqtl$UNIPROT[1]]$note,
    should_create_id = FALSE,
    ID = df_index[consortium == pqtl$study[1] & trait == pqtl$UNIPROT[1]]$id)
}



options(future.globals.maxSize= 1e10)
plan(multisession, workers = 40)

df_index_copy <- df_index
setDT(df_index_copy)
tic()

future_map(split_rosmap, function(x) {brain_pqtl_formatvcf_wrapper( pqtl = x, df_index = df_index_copy )},
           .options = furrr_options(seed = TRUE))

future_map(split_banner, function(x) {brain_pqtl_formatvcf_wrapper( pqtl = x, df_index = df_index_copy )},
           .options = furrr_options(seed = TRUE))

fwrite(df_index, "/mnt/sda/gagelo01/Vcffile/server_gwas_id.txt")

message("this script finished without errors")

toc()

