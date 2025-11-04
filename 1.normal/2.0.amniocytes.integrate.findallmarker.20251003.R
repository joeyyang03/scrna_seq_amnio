rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/1.normal/2.0.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 


library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(ggpubr)
load("/share/home/yangjingyi/project/1.BD/up.load/1.normal/1.0.output/1.0.amniocytes.integrated.harmony.Robj")
p <- DimPlot(object = amniocytes.sum1, reduction = "umap",label = T, repel = TRUE)
##selected cluster
metadata1 <- amniocytes.sum1@meta.data
dim(metadata1)
metadata1$group <- NA
metadata1$group[which(metadata1$seurat_clusters %in% c(0,1,2,3,6))] <- "Cluster1"
metadata1$group[which(metadata1$seurat_clusters %in% c(4,10,14))] <- "Cluster2"
metadata1$group[which(metadata1$seurat_clusters == 5)] <- "Cluster3"
metadata1$group[which(metadata1$seurat_clusters == 7)] <- "Cluster4"
metadata1$group[which(metadata1$seurat_clusters == 8)] <- "Cluster5"
metadata1$group[which(metadata1$seurat_clusters == 9)] <- "Cluster6"
metadata1$group[which(metadata1$seurat_clusters == 11)] <- "Cluster7"
metadata1$group[which(metadata1$seurat_clusters == 12)] <- "Cluster8"
metadata1$group[which(metadata1$seurat_clusters == 13)] <- "Cluster9"
metadata1$group[which(metadata1$seurat_clusters == 15)] <- "Cluster10"


amniocytes.sum1@meta.data <- metadata1
DimPlot(object = amniocytes.sum1, reduction = "umap", group.by = "group", label = T, repel = TRUE)


Idents(amniocytes.sum1) <- "group"
amniocytes.sum1.marker <- FindAllMarkers(amniocytes.sum1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
amniocytes.sum1.marker <- amniocytes.sum1.marker[!duplicated(amniocytes.sum1.marker$gene),]
write.csv(x = amniocytes.sum1.marker, row.names = TRUE, file = "2.0.amniocytes.sum1.slected.cluster.marker.20251003.csv")
