#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggrepel)
setwd("/mnt/sda/gagelo01/Projects/Brain_pQTL/")
dt_res <- fread("Data/Modified/dt_res_pwas.txt")
ntest <- dt_res$UniProt %>% unique %>% length #number of different proteins/different genes tested number of test
res_pwas <- fread("Data/Modified/res_pwas.txt")
susiesign <- fread( file = "Data/Modified/susiesingstringent.txt" )
#Volcano Plot

dt_res <- dt_res[!is.na(b.wald), ]
dt_res[, diffexpressed := "NS"]
dt_res[, diffexpressed := diffexpressed %>% ifelse(b.wald > 0 & pval.wald < 0.05/ntest, "Associated with higher body weight", .) %>%
                          ifelse(b.wald < 0 & pval.wald < 0.05/ntest,  "Associated with lower body weight", .)]
dt_res[, diffexpressed := factor(diffexpressed, levels = c("Associated with higher body weight","Associated with lower body weight", "NS"))]
res_pwas[, delabel := hgnc_symbol]
dt_res <- merge(dt_res, res_pwas[,.(hgnc_symbol, study, delabel)], by = c("hgnc_symbol", "study"), all = TRUE)
dt_res[, logpval := -log10(pval.wald)]


volcano <- ggplot(data = dt_res, aes(x = b.wald, y = logpval, col = diffexpressed,  label=delabel)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values=c("#5884E5","#9E131E", "#B0B0B0")) +
  geom_hline(  yintercept = -log10(0.05/ntest), color = "#A40606", size = 1) +
  geom_label_repel(
    # color = "white",
                   # segment.color = "gray20",
                   # fontface = 'bold',
                   force = 20,
                   max.iter = 1e7,
                   show.legend = F,
                   label.size = 0.5,
                   size = 2,
                   # label.padding = unit(0.2, 'lines'),
                   # direction = "both",
                   # box.padding = 1,
                   # point.padding = unit(0.45, 'lines'),
                   # point.padding = unit(0.5, 'lines'),
                   # min.segment.length = unit(5, 'lines'),
                   # nudge_x = ifelse(subset(merge_MR_res_IVW, pval < bonferroni_threshold)$z_score < 0, (range_y/100*10),(-1*range_y/100*10)),
                   # nudge_y = (range_y/100*10),
                   segment.size = 0.5,
                   max.overlaps = nrow(res_pwas)
  ) +
  labs(x = "Effect on BMI (Beta)", y = expression(-Log[10](P))) +
  theme(    panel.grid.major.y = element_line(size = 0.5, colour = "gray60"),
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
              margin = margin(l = 0.2, r = 0.2)))

ggsave(plot = volcano, filename = "Results/Volcano_pwas_plot.png",
       width = 652/72,height = 585/72,units="in",scale=1, device = "png")


####
dt_res <- fread("Data/Modified/dt_res_pwas.txt")
dt_res <- dt_res[!is.na(beta.multi_cis) & !is.na(pval.multi_cis), ]
dt_res[, diffexpressed := "NS"]
dt_res[, diffexpressed := diffexpressed %>% ifelse(beta.multi_cis > 0 & pval.multi_cis < 0.05/ntest, "Associated with higher body weight", .) %>%
         ifelse(beta.multi_cis < 0 & pval.multi_cis < 0.05/ntest,  "Associated with lower body weight", .)]
dt_res[, diffexpressed := factor(diffexpressed, levels = c("Associated with higher body weight","Associated with lower body weight", "NS"))]

susiesign[, delabel := hgnc_symbol]
dt_res <- merge(dt_res, susiesign[,.(hgnc_symbol, study, delabel)], by = c("hgnc_symbol", "study"), all = TRUE)
dt_res[, logpval := -log10(pval.multi_cis)]
dt_res <- dt_res[beta.multi_cis < 4*sd(beta.multi_cis) & beta.multi_cis > -4*sd(beta.multi_cis),]


volcano <- ggplot(data = dt_res, aes(x = beta.multi_cis, y = logpval, col = diffexpressed,  label=delabel)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values=c("#5884E5","#9E131E", "#B0B0B0")) +
  geom_hline(  yintercept = -log10(0.05/ntest), color = "#A40606", size = 1) +
  geom_label_repel(
    # color = "white",
    # segment.color = "gray20",
    # fontface = 'bold',
    force = 20,
    max.iter = 1e7,
    show.legend = F,
    label.size = 0.5,
    size = 2,
    # label.padding = unit(0.2, 'lines'),
    # direction = "both",
    # box.padding = 1,
    # point.padding = unit(0.45, 'lines'),
    # point.padding = unit(0.5, 'lines'),
    # min.segment.length = unit(5, 'lines'),
    # nudge_x = ifelse(subset(merge_MR_res_IVW, pval < bonferroni_threshold)$z_score < 0, (range_y/100*10),(-1*range_y/100*10)),
    # nudge_y = (range_y/100*10),
    segment.size = 0.5,
    max.overlaps = nrow(susiesign)
  ) +
  labs(x = "Effect on BMI (Beta)", y = expression(-Log[10](P))) +
  theme(    panel.grid.major.y = element_line(size = 0.5, colour = "gray60"),
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
              margin = margin(l = 0.2, r = 0.2)))


ggsave(plot = volcano, filename = "Results/Volcano_pwasmulticis_plot.png",
       width = 652/72,height = 585/72,units="in",scale=1, device = "png")

message("This script finished without errors")
