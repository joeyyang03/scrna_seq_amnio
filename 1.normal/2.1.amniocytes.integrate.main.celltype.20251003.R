rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/1.normal/2.1.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 



library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(ggpubr)
load("/share/home/yangjingyi/project/1.BD/up.load/1.normal/1.0.output/1.0.amniocytes.integrated.harmony.Robj")
amniocytes.sum1.marker <- read.csv("/share/home/yangjingyi/project/1.BD/up.load/1.normal/2.0.output/2.0.amniocytes.sum1.slected.cluster.marker.20251003.csv")
amniocytes.sum1.marker <- amniocytes.sum1.marker%>%as.data.table(.)%>%
  .[(abs(avg_log2FC)>=2)&(pct.1>=0.4|pct.2>=0.4)] #269

table(amniocytes.sum1.marker$cluster)


amniocytes.selected.marker <- amniocytes.sum1.marker %>%
  group_by(cluster) %>% 
  top_n(n = 5, 
        wt = avg_log2FC)#47
amniocytes.selected.marker$gene

  
p1 <- FeaturePlot(amniocytes.sum1, features = c("C1QC","MNDA","HLA-DQA1","MPO",  "HDC", "IL7R"), combine = FALSE)
fix.sc <- scale_color_gradientn(colours = c("#D5D4D4", "#8A0993"),  limits = c(0, 5))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)

p1 <- FeaturePlot(amniocytes.sum1, features = c("KRT78","EPCAM", "COL3A1","COL4A1"), combine = FALSE)
fix.sc <- scale_color_gradientn(colours = c("#D5D4D4", "#2B9070"),  limits = c(0, 5))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)

######################
###marker in paper###
#####################

fetal.liver.hemato.marker <- c("HLA-DRA","LTB", "CD3E","LYZ","CD4", "KIT")#from paper Decoding human fetal liver haematopoiesis
epithelial.cell.marker <- c("EPCAM", "CDH1", "KRT8", "KRT10", "KRT17", "KRT19")#from paper 2024nm for amniotic fluid

marker <- c(fetal.liver.hemato.marker, epithelial.cell.marker)
metadata1 <- amniocytes.sum1@meta.data
metadata1$cell.line <- NA
metadata1$cell.line[which(metadata1$seurat_clusters %in% c(0,1,2,3,6,8,11,15))] <- "Epithelial.cells"
metadata1$cell.line[which(metadata1$seurat_clusters %in% c(4,5,7,9,10,12,13,14))] <- "Hematopoietic.cells"
amniocytes.sum1@meta.data <- metadata1
DimPlot(object = amniocytes.sum1, reduction = "umap", group.by = "cell.line", label = F, repel = TRUE, cols = c("#b2db87","#a48cbe"))

DotPlot(amniocytes.sum1, group.by = 'cell.line',cols = c("white", "#8A0993"),
        features = marker) + RotatedAxis()


amnio.hemato <- subset(x = amniocytes.sum1, subset = cell.line == "Hematopoietic.cells")
table(amnio.hemato@meta.data$cell.line)

save(amnio.hemato, file = "/share/home/yangjingyi/project/1.BD/up.load/1.normal/2.1.output/2.1.amnio.hemato.Robj")

amnio.epi <- subset(x = amniocytes.sum1, subset = cell.line == "Epithelial.cells")
table(amnio.epi@meta.data$cell.line)

save(amnio.epi, file = "/share/home/yangjingyi/project/1.BD/up.load/1.normal/2.1.output/2.1.amnio.epi.Robj")
