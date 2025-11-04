rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/2.1.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 



library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(ggpubr)
load("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/1.1.output/1.1.amniocytes.integrated.harmony.cohort2.Robj")
p <- DimPlot(object = amniocytes.sum1, reduction = "umap",label = T, repel = TRUE)

p1 <- FeaturePlot(amniocytes.sum1, features = c("C1QC","MNDA","HLA-DQA1","MPO",  "HDC", "IL7R"), combine = FALSE )
fix.sc <- scale_color_gradientn(colours = c("#D5D4D4", "#8A0993"),  limits = c(0, 5))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)

p1 <- FeaturePlot(amniocytes.sum1, features = c("KRT78","EPCAM", "COL3A1","COL4A1"), combine = FALSE)
fix.sc <- scale_color_gradientn(colours = c("#D5D4D4", "#2B9070"),  limits = c(0, 5))
p2 <- lapply(p1, function (x) x + fix.sc)
CombinePlots(p2)

###
fetal.liver.hemato.marker <- c("HLA-DRA","LTB", "CD3E","LYZ","CD4", "KIT")
epithelial.cell.marker <- c("EPCAM", "CDH1", "KRT8", "KRT10", "KRT17", "KRT19")

marker <- c(fetal.liver.hemato.marker, epithelial.cell.marker)
metadata1 <- amniocytes.sum1@meta.data
metadata1$cell.line <- NA
metadata1$cell.line[which(metadata1$seurat_clusters %in% c(0,1,2,4,6,7,11,12))] <- "Epithelial.cells"
metadata1$cell.line[which(metadata1$seurat_clusters %in% c(3,5,8,9,10,13))] <- "Hematopoietic.cells"
amniocytes.sum1@meta.data <- metadata1
#####
table(metadata1$cell.line)
 
amniocytes.sum1@meta.data <- metadata1


amnio.hemato <- subset(x = amniocytes.sum1, subset = cell.line == "Hematopoietic.cells")
table(amnio.hemato@meta.data$cell.line)

save(amnio.hemato, file = "/share/home/yangjingyi/project/1.BD/up.load/3.FGR/2.1.output/2.1.cohort2.amnio.hemato.Robj")

amnio.epi <- subset(x = amniocytes.sum1, subset = cell.line == "Epithelial.cells")
table(amnio.epi@meta.data$cell.line)

save(amnio.epi, file = "/share/home/yangjingyi/project/1.BD/up.load/3.FGR/2.1.output/2.1.cohort2.amnio.epi.Robj")

