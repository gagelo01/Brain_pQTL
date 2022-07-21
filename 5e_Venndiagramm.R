#!/usr/bin/env Rscript
#####Venn diagramm
###############Create Venn diagramm #########
library(data.table)
library(tidyverse)
library(VennDiagram)
library(gridExtra)
# your data
setwd("/home/gagelo01/workspace/Projects/Brain_pQTL/")
res_map <- fread("Data/Modified/res_map.txt")
res_pwas <- fread("Data/Modified/res_pwas.txt")
susiesign <- fread("Data/Modified/susiesingstringent.txt")
dt_res <- fread("Data/Modified/dt_res_pwas.txt")

un<- res_map$hgnc_symbol %>% unique
deux<-res_pwas$hgnc_symbol %>% unique
trois<-susiesign$hgnc_symbol %>% unique
# Generate plot
v <- venn.diagram(list('Uni-cis PWMR'=deux, 'Genome-wide mapping'=un, 'Multi-cis PWMR' = trois),
                  fill = c("#F4FAFE", "#4981BF", "#EFC000FF"  ), #"#EFC000FF", "#0073C2FF", "#868686FF", "#CD534CFF"
                  alpha = c(0.5, 0.5, 0.5), main.cex = 4, sub.cex = 0.8, #cex=1.5,
                  main.pos = c(0, 1.05),
                  filename=NULL)


# have a look at the default plot
grid.newpage()
grid.draw(v)

input<-list('Genome-wide \nmapping'=un, 'Uni-cis PWMR'=deux, 'Multi-cis PWMR' = trois)
###########
ven_3_set <- function(input, cex = 1, cat.cex = 1, fontface = "plain", cat.fontface = "plain") {
  func_intern <- function(x) {
    if(length(x) == 0) { return("")}
    if(length(x) == 1) {return(x)}
    if(length(x) > 1) {return(paste(paste(x, c(", ","\n"), sep=""), collapse="") %>%
                                gsub(", $", "", .))}
  }
  # Generate 3 sets of 200 words
  set1 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
  set2 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
  set3 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
  
  # Chart
  v <- venn.diagram(
    x = list(set1, set2, set3),
    category.names = names(input),
    fill = c("#F4FAFE", "#4981BF", "#EFC000FF"  ), #"#EFC000FF", "#0073C2FF", "#868686FF", "#CD534CFF"
    alpha = c(0.5, 0.5, 0.5), main.cex = 12, sub.cex = 12, cex=cex, cat.cex = cat.cex,
    fontface = fontface, cat.fontface = cat.fontface,
    main.pos = c(0, 1.05),
    filename = NULL,
    output=FALSE
  )
  # grid.newpage()
  # grid.draw(v)
  
  # lapply(v, function(i) i$label)
  
  v[[7]]$label <- func_intern(setdiff(input[[1]], unlist(input[2:3])))
  v[[8]]$label <- func_intern(setdiff(intersect(input[[1]], input[[2]]), input[[3]]))
  v[[9]]$label <- func_intern(setdiff(input[[2]], unlist(input[c(1,3)])))
  v[[10]]$label <- func_intern(setdiff(intersect(input[[1]], input[[3]]), input[[2]]))
  v[[11]]$label <- func_intern(intersect(intersect(input[[1]], input[[2]]), input[[3]]))
  v[[12]]$label <- func_intern(setdiff(intersect(input[[2]], input[[3]]), input[[1]]))
  v[[13]]$label <- func_intern(setdiff(input[[3]], unlist(input[1:2])))
  return(v)
}

v <- ven_3_set(input = input, cex = 1.1, cat.cex = 1.2, fontface = "italic", cat.fontface = "bold")
grid.newpage()
grid.draw(v)

# plot  
grid.newpage()
png(filename = "Results/VennDiagramm_comparepwas_mapping.png", width = 680, height = 680,units = "px", pointsize = 12)
grid.draw(v)
# grid.arrange(gTree(children=v), top="Venn Diagramm of protein overlap between Mapping and PWAS")
dev.off()


v <- venn.diagram(list('Ros/Map \n(n = 330)'  = dt_res[study == "rosmap", unique(UniProt)],
                       'Banner \n(n = 149)' = dt_res[study == "banner", unique(UniProt)],
                       'WUSM \n(n = 380)' = dt_res[study == "yang", unique(UniProt)]),
                  fill = c("#F4FAFE", "#4981BF", "#EFC000FF"  ), #"#EFC000FF", "#0073C2FF", "#868686FF", "#CD534CFF"
                  alpha = c(0.5, 0.5, 0.5), main.cex = 4, sub.cex = 4, cex=4, cat.cex = 1,
                  main.pos = c(0, 1.05),
                  filename=NULL)

grid.newpage()
grid.draw(v)
# plot  
grid.newpage()
png(filename = "Results/VennDiagramm_comparestudy.png", width = 437, height = 429,units = "px", pointsize = 16)
grid.draw(v)
# grid.arrange(gTree(children=v), top="Venn Diagramm of protein overlap between Mapping and PWAS")
dev.off()

####
res_all<- rbindlist(list(susiesign, res_pwas, res_map), fill = TRUE)
v <- venn.diagram(list('Ros/Map'  = res_all[study == "rosmap", unique(UniProt)],
                       'Banner' = res_all[study == "banner", unique(UniProt)],
                       'WUSM' = res_all[study == "yang", unique(UniProt)]),
                  fill = c("#F4FAFE", "#4981BF", "#EFC000FF"  ), #"#EFC000FF", "#0073C2FF", "#868686FF", "#CD534CFF"
                  alpha = c(0.5, 0.5, 0.5), main.cex = 4, sub.cex = 4, cex=4, cat.cex = 1,
                  main.pos = c(0, 1.05),
                  filename=NULL)

grid.newpage()
grid.draw(v)
# plot  
grid.newpage()
png(filename = "Results/VennDiagramm_comparestudycausalprot.png", width = 437, height = 429,units = "px", pointsize = 16)
grid.draw(v)
# grid.arrange(gTree(children=v), top="Venn Diagramm of protein overlap between Mapping and PWAS")
dev.off()

message("This script finished without errors")
