#!/usr/bin/env Rscript
library(GagnonMR)
library(biomaRt)
library(data.table)
library(tidyverse)
library(CMplot)
library(ggrepel)
library(gassocplot)
library(cowplot)
library(ggplot2)
library(viridis)
library(VennDiagram)
library(ggpubr)

setwd("/home/gagelo01/workspace/Projects/Brain_pQTL")
gencode = data.table::fread("/home/couchr02/workspace/GTEx_v8/gencode.v26.GRCh38.genes.txt", data.table = F, stringsAsFactors = F)
res_map <- fread( "Data/Modified/res_map.txt")
res_pwas <- fread( "Data/Modified/res_pwas.txt")
dt_gene_region <- fread( "Data/Modified/dt_gene_region.txt")
tpm <- fread( "Data/Modified/data_united_tpm.txt")
dt_res <- fread("Data/Modified/dt_res_pwas.txt")
ntest <- dt_res[!is.na(b.wald), ]$UniProt %>% unique %>% length #number of different proteins tested number of test
res_hypr <- fread( "Data/Modified/hypr_res.txt")
exposure_hypr <- fread( "Data/Modified/exposure_hypr.txt")
res_specificity <- fread("Data/Modified/res_specificity.txt")

#######Fig 1
source("Analysis/CMplot_eloi.R")
gwasvcf::set_bcftools()
gwasvcf::set_plink()
tomanhatthan<- res_map[ , .SD[which.max(posprob_coloc_PPH4)] , by = "hgnc_symbol", .SDcols = c("posprob_colocH4.SNP")]
data <- gwasvcf::query_gwas("/mnt/sda/gagelo01/Vcffile/Server_vcf/trait-1-1/trait-1-1.vcf.gz", pval = 5e-5) %>%
  gwasglue::gwasvcf_to_TwoSampleMR(., type = "outcome") %>% as.data.table(.)

data <- data[!(is.na(chr.outcome) | is.na(pos.outcome) | is.na(pval.outcome)) & pval.outcome < 5e-5,]
data <- data[,.(SNP, chr.outcome, pos.outcome, pval.outcome)]

setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL/Results")

CMplot_eloi(as.data.frame(data), type="p", plot.type="m", LOG10=TRUE, threshold=5e-08, amplify=FALSE, memo="", dpi=300, verbose=TRUE, width=14,
            height=6, col=c("grey70", "grey90"), threshold.lwd=2, cex=0.4, highlight.cex =  0.6, highlight=tomanhatthan$posprob_colocH4.SNP,
            highlight.text=tomanhatthan$hgnc_symbol, highlight.col = "red",
            highlight.text.xadj =rep(-1, length(tomanhatthan$hgnc_symbol)), highlight.text.yadj =rep(1, length(tomanhatthan$hgnc_symbol)),
            file.output = TRUE, arrow_col_eloi = "black", highlight.text.cex=1, padding = 3)

setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL/")
#######Volcano#################
volcano <- dt_res[!is.na(b.wald), ]
volcano[, diffexpressed := "NS"]
volcano[, diffexpressed := diffexpressed %>% ifelse(b.wald > 0 & pval.wald < 0.05/ntest, "Associated with higher body weight", .) %>%
         ifelse(b.wald < 0 & pval.wald < 0.05/ntest,  "Associated with lower body weight", .)]
volcano[, diffexpressed := factor(diffexpressed, levels = c("Associated with higher body weight","Associated with lower body weight", "NS"))]
res_pwas[pval.wald < 0.05/ntest & posprob_coloc_PPH4 > 0.8  & res_pwas$steiger_dir.wald == TRUE & fstat.exposure.wald > 10, delabel := hgnc_symbol]
volcano <- merge(volcano, res_pwas[,.(id.exposure, delabel)], by = "id.exposure", all.x = TRUE)
volcano[, logpval := -log10(pval.wald)]
volcano <- volcano[steiger_dir.wald==TRUE,]
volcano <- volcano[!(hgnc_symbol == "EIF5B" & study == "rosmap"),]
volcano[is.na(delabel), diffexpressed := "Not prioritised"]


lim <- ifelse(test = abs(min(volcano$b.wald)) > abs(max(volcano$b.wald)),
              yes = abs(min(volcano$b.wald)),
              no = abs(max(volcano$b.wald)))

plot_volcano <- function(volcano, xlab = "Effect on BMI (Beta)", lim = 1) {
  
  # lim_y1 <- mean(volcano$logpval) * 5 * sd(volcano$logpval)
  # lim_y2 <- max(volcano[!is.na(delabel),]$logpval)
  # lim_y <- max(c(lim_y1, lim_y2))
  
  lim_y <- max(volcano$logpval)
  volcanoplot <-  ggplot(data = volcano, aes(x = b.wald, y = logpval, col = diffexpressed,  label=delabel)) +
    geom_point(size = 1, alpha = ifelse(volcano[, is.na(delabel)], 0.5, 1)) +
    theme_bw() +
    scale_color_manual(values=c("#5884E5","#9E131E", "#B0B0B0")) +
    geom_hline(  yintercept = -log10(0.05/ntest), color = "#A40606", size = 1) +
    ggrepel::geom_text_repel(
      force = 5,
      force_pull = 1,
      box.padding = 0.4,
      max.iter = 1e7,
      show.legend = F,
      text.size = 1.5,
      linewidth = 3.5,
      segment.size = 0.5,
      segment.alpha = 0.5,
      segment.linetype = "solid",
      min.segment.length = 0,
      max.overlaps = nrow(volcano[!is.na(delabel), ]) +10,
      nudge_y = volcano$logpval,
      nudge_x = volcano$b.wald*2#(abs(volcano$b.wald)+lim/10) * volcano$b.wald/abs(volcano$b.wald)
    ) +
    labs(x = xlab, y = expression(-Log[10](P))) +
    theme(    panel.grid.major.y = element_line(linewidth = 0.5, colour = "gray60"),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.background = element_blank(),
              plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, "cm"),
              legend.position = "top",
              legend.title = element_blank(),
              # axis.title = element_text(size = 25, colour = "gray20"),
              axis.line = element_line(size = 1, colour = "gray20"),
              axis.ticks = element_line(size = 1, colour = "gray20"),
              axis.ticks.length = unit(.25, "cm"),
              legend.text = element_text(
                color = "gray20",
                size = 10,
                margin = margin(l = 0.2, r = 0.2))) +
    scale_x_continuous(limits=c(-lim, lim)) +
    scale_y_continuous(limits=c(0, lim_y)) +
    geom_vline(xintercept = 0, lty = 2)
  
  return(volcanoplot)
}

volcanoplot <- plot_volcano(volcano = volcano, lim = lim)
volcanoplot
ggsave(plot = volcanoplot, filename = "Results/Volcano_pwas_plot.png",
       width = 652/72,height = 585/72,units="in",scale=1, device = "png")

# volcanoplot <- ggplot(data = volcano, aes(x = b.wald, y = logpval, col = diffexpressed,  label=delabel)) +
#   geom_point() +
#   theme_bw() +
#   scale_color_manual(values=c("#5884E5","#9E131E", "#B0B0B0")) +
#   geom_hline(  yintercept = -log10(0.05/ntest), color = "#A40606", size = 1) +
#   geom_label_repel(
#     # color = "white",
#     # segment.color = "gray20",
#     # fontface = 'bold',
#     force = 20,
#     max.iter = 1e7,
#     show.legend = F,
#     label.size = 0.5,
#     size = 2,
#     # label.padding = unit(0.2, 'lines'),
#     # direction = "both",
#     # box.padding = 1,
#     # point.padding = unit(0.45, 'lines'),
#     # point.padding = unit(0.5, 'lines'),
#     # min.segment.length = unit(5, 'lines'),
#     # nudge_x = ifelse(subset(merge_MR_res_IVW, pval < bonferroni_threshold)$z_score < 0, (range_y/100*10),(-1*range_y/100*10)),
#     # nudge_y = (range_y/100*10),
#     segment.size = 0.5,
#     max.overlaps = nrow(res_pwas)
#   ) +
#   labs(x = "Effect on BMI (Beta)", y = expression(-Log[10](P))) +
#   theme(    panel.grid.major.y = element_line(size = 0.5, colour = "gray60"),
#             panel.grid.major.x = element_blank(),
#             panel.grid.minor.y = element_blank(),
#             panel.grid.minor.x = element_blank(),
#             panel.background = element_blank(),
#             plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, "cm"),
#             legend.position = "top",
#             legend.title = element_blank(),
#             # axis.title = element_text(size = 25, colour = "gray20"),
#             axis.line = element_line(size = 1, colour = "gray20"),
#             axis.ticks = element_line(size = 1, colour = "gray20"),
#             axis.ticks.length = unit(.25, "cm"),
#             legend.text = element_text(
#               color = "gray20",
#               size = 10,
#               margin = margin(l = 0.2, r = 0.2)))
# 
# volcanoplot
# 
# ggsave(plot = volcanoplot, filename = "Results/Volcano_pwas_plot.png",
#        width = 652/72,height = 585/72,units="in",scale=1, device = "png")

###Venn diagramm
un<- res_map$hgnc_symbol %>% unique
deux<-res_pwas$hgnc_symbol %>% unique
input<-list('Genome-wide \nmapping'=un, 'Uni-cis PWMR'=deux)

ven_2_set <- function(input, cex = 1, cat.cex = 1, fontface = "plain", cat.fontface = "plain", word_per_line = 3) {
  
  func_intern <- function(x, word_per_line) {
    if(length(x) == 0) { return("")}
    if(length(x) == 1) {return(x)}
    if(length(x) > 1) {return(paste(paste(x, c(rep("  ", word_per_line -1),"\n"), sep=""), collapse="") %>%
                                gsub(", $", "", .))}
  }
  # Generate 3 sets of 200 words
  set1 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
  set2 <- paste(rep("word_" , 100) , sample(c(1:1000) , 100 , replace=F) , sep="")
  
  # Chart
  v <- venn.diagram(
    x = list(set1, set2),
    category.names = names(input),
    fill = c("#4981BF", "#EFC000FF"), #"#EFC000FF", "#0073C2FF", "#868686FF", "#CD534CFF"
    alpha = c(0.5), main.cex = 16, sub.cex = 18, cex=cex, cat.cex = cat.cex,
    fontface = fontface, cat.fontface = cat.fontface, main.fontfamily = "Helvetica Sans",
    sub.fontfamily = "Helvetica Sans", fontfamily = "Helvetica Sans",
    cat.pos = c(330, 30),
    main.pos = c(0, 1.05),
    filename = NULL,
    output=FALSE,
    scaled = FALSE
  )
  # grid.newpage()
  # grid.draw(v)
  
  # lapply(v, function(i) i$label)
  
  v[[5]]$label <- func_intern(setdiff(input[[1]], input[[2]]), word_per_line = word_per_line[1])
  v[[6]]$label <- func_intern(setdiff(input[[2]], input[[1]]), word_per_line = word_per_line[3])
  v[[7]]$label <- func_intern(intersect(input[[2]], input[[1]]), word_per_line = word_per_line[2])
  
  return(v)
}

v <- ven_2_set(input = input, cex = 0.80, cat.cex = 0.95,
               fontface = "italic", cat.fontface = "bold", word_per_line = c(2,2,3))
grid.newpage()
grid.draw(v)
# plot  
reso <- 1200
length <- 480*reso/72
grid.newpage()
png(file="Results/VennDiagramm_comparepwas_mapping.png", units="px", 
    res=reso,height=length,width=length ) # or other device
# png(filename = "Results/VennDiagramm_comparepwas_mapping.png", width = 500, height = 500,units = "px", pointsize = 12)
grid.draw(v)
# grid.arrange(gTree(children=v), top="Venn Diagramm of protein overlap between Mapping and PWAS")
dev.off()


#####hyprcoloc plot
dir.create("Results/Hypr_plot")
exposure_vec <- c("ADCY3", "DOC2A")
for(i in 1:length(exposure_vec)) {
  
  A <- stack_assoc_plot_wrapper(df_aligned = exposure_hypr[gene.exposure == exposure_vec[i]],
                                res_hypr1 = res_hypr[gene_investigated== exposure_vec[i]],
                                ldref = "/home/couchr02/Mendel_Commun/Christian/LDlocal/EUR_rs",
                                traits_inorder = unique(exposure_hypr[gene.exposure == exposure_vec[i]]$exposure),
                                build = 37)
  res <- sensitivity.plot_wrapper(df_aligned = exposure_hypr[gene.exposure == exposure_vec[i]],
                                  traits_inorder = unique(exposure_hypr[gene.exposure == exposure_vec[i]]$exposure))
  B<-drawheatmap(res[[2]])
  twopanel <-  cowplot::ggdraw() +
    draw_plot(ggplotify::as.ggplot(A) + theme(text = element_text(size = 0.4)), x = 0.08, y =0, width = .6, height = 1) +
    draw_plot(B, x = .66, y =0.1, width = .31, height = 0.7) +
    draw_plot_label(label = c("", ""), size = 25,
                    x = c(0, 0.62), y = c(0.9, 0.9))
saveRDS(object = twopanel, file = paste0("Results/Hypr_plot/", "twopanel_hypr_plot_",exposure_vec[i], ".rds"))
  ggsave(plot = twopanel, filename = paste0("Results/Hypr_plot/", "twopanel_hypr_plot_",exposure_vec[i], ".png"),
         width = 790/72,height = 683/72,units="in",scale=1, device = "png")
}


list_plot <- map(exposure_vec, function(x) readRDS(file = paste0("Results/Hypr_plot/", "twopanel_hypr_plot_",x, ".rds")))

ggpubr::ggarrange(list_plot[[1]], list_plot[[2]], 
          labels = c("A)", "B)"),
          ncol = 2, nrow = 1)

ggsave(filename =  paste0("Results/twopanel_hypr_plot_AB.png"),  
       width = 1535/72, height = 775/72, units = "in", scale = 1, device="png" )
#######histogramm tissue celltype############
# histogramm_tissue_celltype <- function(M) {
#   M[,logp := -log10(P)]
#   M[, tissue_name := gsub("_", " ", FULL_NAME)]
#   M <- M[order(-logp)]
#   M[,tissue_name := factor(tissue_name, levels = unique(tissue_name))]
#   bonfp <- M[,-log10(0.05/.N)] #bonferroni correction threshold
#   nomp <- -log10(0.05)
#   M[, colourbar := logp %>% ifelse(. > bonfp, "red3", .) %>% ifelse(logp > nomp & logp < bonfp, "blue", .) %>% ifelse(. != "red3" & . != "blue", "grey", .)]
#   M[,colourbar := factor(colourbar, levels = c("red3", "blue", "grey"))]
#   
#   ggplot(M, aes(x = tissue_name, y =logp, fill = colourbar)) +
#     geom_bar(stat = "identity") +
#     theme(axis.text.x = element_text(face = "bold", angle = 85, hjust = 1)) +
#     scale_fill_manual("legend", values = levels(M$colourbar)[levels(M$colourbar) %in% unique(M$colourbar)]) +
#     theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
#     geom_hline(yintercept = c(bonfp, nomp), linetype = "dashed", color = c("red3", "blue")) + 
#     guides(fill="none") +
#     ylab("-log 10 P-Value") +
#     xlab("")
#   
# }
# 
# histogramm_tissue_celltype(M = tissue)
# 
# ggsave("Results/barplot_tissular_specificiy_.png", 
#        width = 831/72,
#        height = 466/72,
#        units = "in")
# 
# #Cell type
# 
# histogramm_tissue_celltype(M = celltype)
# 
# ggsave("Results/barplot_celltype_enrichment.png",
#        width=290/72,
#        height=250/72,
#        units = "in")

#UpSet plot to evaluate gene overlap
res_specificity <- res_specificity[UniProt%in%res_pwas$UniProt,]
k <- rbindlist(list(res_pwas, res_specificity), fill = TRUE)
listInput <- lapply(split(k, k$study), function(x) unique(x$UniProt))
plotupset <- UpSetR::upset(UpSetR::fromList(listInput), order.by = "freq", nsets = length(listInput))
png(file="Results/upset_full.png") # or other device
plotupset
dev.off()

kcause <- k[pval.wald < 0.05/ntest & posprob_coloc_PPH4 > 0.8  & steiger_pval.wald < 0.05 & fstat.exposure.wald > 10, unique(UniProt), by = "study"]
studyvec<-c(yang="Yang et al.", deCODE = "deCODE", FENLAND = "Fenland", GTEX = "GTEx", eQTLGen = "eQTLGen",
            banner = "Banner", rosmap = "ROS/MAP")
kcause[, study := studyvec[study]]
listInput <- lapply(split(kcause, kcause$study), function(x) x$V1)
plotupset_cause<-UpSetR::upset(UpSetR::fromList(listInput), order.by = "freq", nsets = length(listInput))
reso <- 1200
length <- 480*reso/72
png(file="Results/upset_cause.png", units="px", res=reso,height=length,width=length ) # or other device
plotupset_cause
dev.off()
######heatmap tissular specificity 
data_united<-tpm
heatmap_tissular_specificity <- function(data_united ) {
  # sorted_gene = unique(dplyr::arrange(data_united, chromosome_name, band) %>% dplyr::select(Description, chr_band))
  # x = rev(sorted_gene$Description)[1]
  label_names = sapply(data_united$Description, FUN = function(x){bquote(paste(.(x), " (",tau,": ", .(format(unique(data_united$fTau[which(data_united$Description == x)]),  digits = 2, nsmall = 2)), ")"))})



  g <- ggplot(data_united, aes(x = tissue, y = factor(Description, levels = rev(unique(data_united$Description))), fill = TPM))  +
    geom_tile() +
    labs(fill = "Averaged TPM\n") +
    scale_fill_viridis(discrete = FALSE) +
    scale_y_discrete(labels = label_names) +
    # guides(colour = guide_colourbar(title.vjust = 1)) +
    theme(
      # legend.position = "none",
      # panel.border = element_rect(colour = "gray20", fill = "transparent", size = 1),
      # panel.grid.major.y = element_line(size = 0.5, colour = "gray60"),
      # panel.grid.major.x = element_blank(),
      # panel.grid.minor.y = element_blank(),
      # panel.grid.minor.x = element_blank(),
      panel.background = element_blank(),
      # panel.grid.major = element_blank(),
      plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, "cm"),
      #legend.position = "right",
      legend.position = "top",
      legend.title = element_text(
        color = "gray20",
        size = 12
        # margin = margin(l = 0.2, r = 0.2, t = 0.2, b = 0.8)
      ),
      legend.text = element_text(
        color = "gray20",
        size = 10
        # margin = margin(l = 0.2, r = 0.2)
      ),
      legend.title.align = 0.5,
      legend.spacing.y = unit(0.1, 'cm'),
      legend.key = element_rect(fill = "transparent", colour = "transparent"),
      legend.key.size = unit(0.8, "cm"),
      # panel.grid.major.x = element_blank(),
      axis.title = element_blank(),
      axis.line = element_line(size = 1, colour = "gray20"),
      axis.ticks = element_line(size = 1, colour = "gray20"),
      axis.text.y = element_text(
        # angle = 60,
        size = 10,
        # vjust = 0.5,
        colour = "gray20"
      ),
      axis.text.x = element_text(
        angle = 60,
        size = 10,
        # vjust = 0.55,
        hjust = 1,
        colour = "gray20"
      ),
      # axis.text = element_text(
      #   # angle = 60,
      #   size = 20,
      #   # vjust = 0.5,
      #   colour = "gray20"
      # ),
      axis.ticks.length = unit(.25, "cm"))

  print(g)
}

setnames(data_united, "Averaged.TPM", "TPM")
list_gene <- list(gene_pwas = setdiff(res_pwas$hgnc_symbol, res_map$hgnc_symbol) %>% unique(),
                  gene_map = setdiff(res_map$hgnc_symbol, res_pwas$hgnc_symbol) %>% unique(),
                  intersection = res_map[hgnc_symbol %in% res_pwas$hgnc_symbol, unique(hgnc_symbol)],
                  all_gene = data_united[order(-fTau),][,unique(Description)])

for(i in 1:length(list_gene)) {
  heatmap_tissular_specificity(data_united[order(-fTau),][Description %in% list_gene[[i]],]) +
    theme(axis.text.x= element_text(face=ifelse(grepl("brain", tolower(unique(data_united$tissue))),"bold","italic")))
  ggsave(paste0("Results/heatmap_tissular_specificiy_", names(list_gene)[i] ,".png"),
         width = 978/72,
         height = (50+length(list_gene[[i]])*25)/72,
         units = "in",
         type = "cairo")
}

message("This script finished without errors")


