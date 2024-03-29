#!/usr/bin/env Rscript
library(GagnonMR)
library(dplyr)
library(biomaRt)
library(data.table)
library(tidyverse)

setwd("/home/gagelo01/workspace/Projects/Brain_pQTL")
gencode = data.table::fread("/home/couchr02/workspace/GTEx_v8/gencode.v26.GRCh38.genes.txt", data.table = F, stringsAsFactors = F)
res_map <- fread( "Data/Modified/res_map.txt")
res_pwas <- fread( "Data/Modified/res_pwas.txt")
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")

get_tpm_for_genes_on_all_tissues <- function(gene_list, 
                                             gencode = data.table::fread("/home/couchr02/workspace/GTEx_v8/gencode.v26.GRCh38.genes.txt", data.table = F, stringsAsFactors = F, nThread = 6),
                                             tissue_gtex_wd =  "/home/couchr02/workspace/GTEx_v8/list/",
                                             gtex_tpm_file_name = "/home/couchr02/Mendel_Commun/GTEx_v8/eQTL/exp_TPM/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct",
                                             individus_EUR = data.table::fread("/home/couchr02/Mendel_Commun/Nicolas/GTEx_V8/EUR_GTEx8.txt", stringsAsFactors = F, data.table = F, header = F)
){
  if (exists("data_tau")){
    rm(data_tau)
  }
  
  if(gene_list[1]!="Name"){gene_list <- c("Name",gene_list)}
  gene_loop = unique(gencode$gene_id[which(gencode$gene_name %in% gene_list[c(2:length(gene_list))])])
  
  data.table::fwrite(list(gene_list), "grep_list.txt", quote = F, col.names = F, row.names = F, sep = "\t")
  
  tpm_count = data.table::fread(cmd = paste0("zgrep -w -f grep_list.txt ", gtex_tpm_file_name),
                                data.table = F,
                                stringsAsFactors = F)
  tpm_count = tpm_count[which(tpm_count$Description %in% gene_list[-1]),]   # gene_list[-1] pour enlever "Name" de la liste
  
  if (file.exists("grep_list.txt")){
    file.remove("grep_list.txt")
  }
  
  for (gene in gene_loop) {
    gene_name = unique(gencode$gene_name[which(gencode$gene_id == gene)])
    message("\n#### Gene in progress: ", gene_name," (", which(gene_loop == gene), "/", length(gene_loop), ") ####")
    
    if (length(which(tpm_count$Description == gene_name)) > 0){
      tpm_count_gene = subset(tpm_count, Description == gene_name)  
    } else {
      message("   Missing gene")
      next()
    }
    
    gene_id = tpm_count_gene$Name
    
    if (gene == gene_loop[1]){
      data_tau = data.frame(Ensembl.Gene.ID = gene, stringsAsFactors = F)
    } else {
      data_tau = rbind(data_tau, NA)
      data_tau$Ensembl.Gene.ID[which(is.na(data_tau$Ensembl.Gene.ID) == T)] = gene
    }
    
    
    
    tissues_loop = list.files(tissue_gtex_wd, pattern="_liste.txt")  
    tissues_loop = tissues_loop[-which(tissues_loop %in% c("Vagina_liste.txt", "Uterus_liste.txt", "Testis_liste.txt", "Prostate_liste.txt", "Ovary_liste.txt", "Breast_Mammary_Tissue_liste.txt", "Fallopian_Tube_liste.txt", "Cervix_Ectocervix_liste.txt", "Cervix_Endocervix_liste.txt", "Bladder_liste.txt", "Kidney_Cortex_liste.txt", "Kidney_Medulla_liste.txt"))]
    
    #tissue_list = tissues_loop[1]
    for (tissue_list in tissues_loop) {
      tissue = gsub(tissue_list, pattern = "_liste.txt", replacement = "")
      
      # message("     ## Tissue in progress: ", tissue," (", which(tissues_loop == tissue_list), "/", length(tissues_loop), ") ####")
      
      list_tissue = data.table::fread(paste0(tissue_gtex_wd, tissue_list), data.table = F, stringsAsFactors = F, nThread = 6, header = F)
      
      list_tissue$ID = sapply(list_tissue$V1, FUN = function(x) {paste(unlist(strsplit(x, "-"))[1], unlist(strsplit(x, "-"))[2], sep = "-")})
      
      list_tissue_EUR = subset(list_tissue, ID %in% individus_EUR$V1)
      
      tpm_count_tissue_EUR = tpm_count_gene[,c("Name", "Description", colnames(tpm_count_gene)[which(colnames(tpm_count_gene) %in% list_tissue_EUR$V1)])]
      
      tpm_count_tissue_EUR = t(tpm_count_tissue_EUR)
      tpm_count_tissue_EUR = tpm_count_tissue_EUR[grep("GTEX-", rownames(tpm_count_tissue_EUR)),, drop = F]
      
      if (gene == gene_loop[1]) {
        data_tau = cbind(data_tau, data.frame(tissue = mean(as.numeric(tpm_count_tissue_EUR[,1])), stringsAsFactors = F))
        colnames(data_tau)[which(colnames(data_tau) == "tissue")] = paste0("Averaged.TPM.", tissue)
      } else {
        data_tau[data_tau$Ensembl.Gene.ID == gene_id, paste0("Averaged.TPM.", tissue)] = mean(as.numeric(tpm_count_tissue_EUR[,1]))
      }
      
    } # end tissues loop
  }# end genes loop
  
  
  tpm = 1
  x = data_tau[,c(-1)]
  x[x < tpm] = 1
  data_tau[,c(-1)] = log2(x)
  
  fTau_TPM <- function(average_TPM)
  {
    if (all(!is.na(average_TPM)))
    {
      if (min(average_TPM, na.rm = TRUE) >= 0)
      {
        if (max(average_TPM) != 0)
        {
          x =  1-(average_TPM/max(average_TPM))
          res <- sum(x, na.rm = TRUE)
          res <- res/(length(x) - 1)
        } else {
          res <- 0
        }
      } else {
        res <- NA
        print("Expression values have to be positive!")
      } 
    } else {
      res <- NA
      print("No data for this gene avalable.")
    } 
    return(res)
  }
  
  
  data_tau$fTau = apply(data_tau[,grep("Averaged.TPM.", colnames(data_tau))], 1, fTau_TPM)
  data_tau = merge(data_tau, tpm_count[,c("Description", "Name")], by.x = "Ensembl.Gene.ID", by.y = "Name")
  
  #### Heatmap ####
  
  # listMarts()
  #ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org",dataset = "hsapiens_gene_ensembl")
  ensembl = useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", version = "91" )
  
  # biomaRt::listEnsembl()
  # listDatasets(useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org"))
  
  filters = listFilters(ensembl)
  
  attributes = listAttributes(ensembl)
  
  # attributes[grep("band", attributes$name),]
  
  biomart_output =  getBM(
    attributes = c("ensembl_gene_id_version", "band", "chromosome_name"),
    filters = "ensembl_gene_id_version",
    values = unique(data_tau$Ensembl.Gene.ID),
    mart = ensembl
  )
  
  data_tau = merge(data_tau, biomart_output, by.x = "Ensembl.Gene.ID", by.y = "ensembl_gene_id_version")
  
  data_tau$chr_band = paste0(data_tau$chromosome_name, data_tau$band)
  
  
  for (i in grep("tpm", colnames(data_tau), ignore.case = T)){
    message("i= ", i)
    if (i == grep("tpm", colnames(data_tau), ignore.case = T)[1])
    {
      data_united = cbind(data_tau[,c(which(colnames(data_tau) %in% c("Description", "band", "chromosome_name", "chr_band", "fTau")), i)], tissue = rep(gsub(colnames(data_tau)[i], pattern = "Averaged.TPM.", replacement = ""), dim(data_tau)[1]))
      colnames(data_united)[which(colnames(data_united) == colnames(data_tau)[i])] = "Averaged.TPM"
    } else {
      data_to_bind = cbind(data_tau[,c(which(colnames(data_tau) %in% c("Description", "band", "chromosome_name", "chr_band", "fTau")), i)], tissue = rep(gsub(colnames(data_tau)[i], pattern = "Averaged.TPM.", replacement = ""), dim(data_tau)[1]))
      colnames(data_to_bind)[which(colnames(data_to_bind) == colnames(data_tau)[i])] = "Averaged.TPM"
      data_united = rbind(data_united, data_to_bind)
    }
  }
  
  data_united$tissue = gsub(data_united$tissue, pattern = "_", replacement = " ")
  data_united$tissue = gsub(data_united$tissue, pattern = "([a-z])([A-Z])", replacement = "\\1 \\2", perl = T)
  
  return(data_united)
}

data_united <- get_tpm_for_genes_on_all_tissues(gene_list = dt_gene_region[UniProt %in% (c(res_map$UniProt, res_pwas$UniProt) %>% unique), hgnc],
                                                gencode = gencode)

fwrite(data_united, "Data/Modified/data_united_tpm.txt")


#####
tpm <- fread( "Data/Modified/data_united_tpm.txt")
data_united <- tpm
#####
toselect <- res_map[, unique(hgnc_symbol)]
setnames(tpm, "Averaged.TPM", "TPM")
tpm[, Averaged.TPM := mean(TPM), by = "Description"]
tpm[, Relative.TPM := TPM/Averaged.TPM]
tpm[ grepl("brain", tolower(tissue)), mean(Relative.TPM), by = "Description"][Description %in% toselect][V1 > 1,][order(-V1)]
###very simply is mean tpm in brain higher than all other tissue combined

brain<- tpm[grepl("brain", tolower(tissue)),  mean(TPM), by = Description] 
setnames(brain, "V1", "brain")
other<- tpm[!grepl("brain", tolower(tissue)), mean(TPM), by = Description]
setnames(other, "V1", "other")

compare<-merge(brain, other, by = "Description")
compare <- compare[, relative := brain/other][order(-relative)]
fwrite(compare, "Data/Modified/compare_tpm.txt")
fwrite(data.frame(Gene = compare$Description), file = "Data/Modified/inputmetascape.csv", sep = ",")
######heatmap
# library(ggplot2)
# library(viridis)
# heatmap_tissular_specificity <- function(data_united ) {
#   # sorted_gene = unique(dplyr::arrange(data_united, chromosome_name, band) %>% dplyr::select(Description, chr_band))
#   # x = rev(sorted_gene$Description)[1]
#   label_names = sapply(data_united$Description, FUN = function(x){bquote(paste(.(x), " (",tau,": ", .(format(unique(data_united$fTau[which(data_united$Description == x)]),  digits = 2, nsmall = 2)), ")"))})
#   
#   
#   
#   g <- ggplot(data_united, aes(x = tissue, y = factor(Description, levels = rev(unique(data_united$Description))), fill = TPM))  +
#     geom_tile() +
#     labs(fill = "Averaged TPM\n") +
#     scale_fill_viridis(discrete = FALSE) +
#     scale_y_discrete(labels = label_names) +
#     # guides(colour = guide_colourbar(title.vjust = 1)) +
#     theme(
#       # legend.position = "none",
#       # panel.border = element_rect(colour = "gray20", fill = "transparent", size = 1),
#       # panel.grid.major.y = element_line(size = 0.5, colour = "gray60"),
#       # panel.grid.major.x = element_blank(),
#       # panel.grid.minor.y = element_blank(),
#       # panel.grid.minor.x = element_blank(),
#       panel.background = element_blank(),
#       # panel.grid.major = element_blank(),
#       plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, "cm"),
#       #legend.position = "right",
#       legend.position = "top",
#       legend.title = element_text(
#         color = "gray20",
#         size = 12
#         # margin = margin(l = 0.2, r = 0.2, t = 0.2, b = 0.8)
#       ),
#       legend.text = element_text(
#         color = "gray20",
#         size = 10
#         # margin = margin(l = 0.2, r = 0.2)
#       ),
#       legend.title.align = 0.5,
#       legend.spacing.y = unit(0.1, 'cm'),
#       legend.key = element_rect(fill = "transparent", colour = "transparent"),
#       legend.key.size = unit(0.8, "cm"),
#       # panel.grid.major.x = element_blank(),
#       axis.title = element_blank(),
#       axis.line = element_line(size = 1, colour = "gray20"),
#       axis.ticks = element_line(size = 1, colour = "gray20"),
#       axis.text.y = element_text(
#         # angle = 60,
#         size = 10,
#         # vjust = 0.5,
#         colour = "gray20"
#       ),
#       axis.text.x = element_text(
#         angle = 60,
#         size = 10,
#         # vjust = 0.55,
#         hjust = 1,
#         colour = "gray20"
#       ),
#       # axis.text = element_text(
#       #   # angle = 60,
#       #   size = 20,
#       #   # vjust = 0.5,
#       #   colour = "gray20"
#       # ),
#       axis.ticks.length = unit(.25, "cm"))
#   
#   print(g)
# }
# 
# list_gene <- list(gene_pwas = setdiff(res_pwas$hgnc_symbol, res_map$hgnc_symbol) %>% unique(),
# gene_map = setdiff(res_map$hgnc_symbol, res_pwas$hgnc_symbol) %>% unique(),
# intersection = res_map[hgnc_symbol %in% res_pwas$hgnc_symbol, unique(hgnc_symbol)],
# susie_list = unique(susiesign$hgnc_symbol),
# all_gene = data_united[order(-fTau),][,unique(Description)])
# 
# for(i in 1:length(list_gene)) {
# heatmap_tissular_specificity(data_united[order(-fTau),][Description %in% list_gene[[i]],]) +
#   theme(axis.text.x= element_text(face=ifelse(grepl("brain", tolower(unique(data_united$tissue))),"bold","italic")))
# ggsave(paste0("Results/heatmap_tissular_specificiy_", names(list_gene)[i] ,".png"), 
#        width = 978/72,
#        height = (length(list_gene[[i]])*25)/72,
#        units = "in",
#        type = "cairo")
# }
# 
# message("This script finished without errors")
# 
# 
