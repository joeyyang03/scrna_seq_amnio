rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/4.1.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 

library(Seurat)
library(ggplot2)
library(tidyverse)
library(rhdf5)
library("edgeR")
library(dplyr)
library(readr)
library(harmony)
library(SingleR)
library(SingleCellExperiment)
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]


cellline <- basename(input_file)
cellline <- gsub(".integrated.20251007.RData", "", cellline)
cellline <- gsub("4.1.cohort2.", "", cellline)

load(input_file)

DefaultAssay(combined_seurat) <- "integrated"
combined_seurat <- ScaleData(object = combined_seurat, verbose = FALSE)
combined_seurat <- RunPCA(object = combined_seurat, npcs = 13, verbose = FALSE)
combined_seurat <- RunUMAP(object = combined_seurat, reduction = "pca", dims = 1:13)
combined_seurat <- FindNeighbors(object = combined_seurat, reduction = "pca", dims = 1:13)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)

ref <- subset(combined_seurat, subset = group %in% c("amniocytes"))
amnio <- subset(combined_seurat, subset = group %!in% c("amniocytes"))



amnio$seurat_clusters <- as.factor(amnio$seurat_clusters)
levels(amnio$seurat_clusters)

results <- SingleR(
  test = amnio@assays$integrated@data, 
  ref = ref@assays$integrated@data,  
  labels = ref$singleR_labels,  
  clusters = amnio$seurat_clusters  
)

###
amnio$singleR_labels <- results@listData$labels[amnio$seurat_clusters]

########
save.image(file = paste0("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/4.1.output/4.1.cohort2.",cellline,".integrated.singleR.with.inte.cluster.res.0.5.20251007.RData"))


#####
rm(list = ls())
load("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/4.1.output/4.1.cohort2.epi.integrated.singleR.with.inte.cluster.res.0.5.20251007.RData")
amnio.epi.metadata <- amnio@meta.data
levels(amnio.epi.metadata$seurat_clusters)
levels(amnio.epi.metadata$integrated_snn_res.0.5)
table(amnio.epi.metadata$RNA_snn_res.0.4)
amnio.epi.metadata$barcode <- rownames(amnio.epi.metadata)
amnio.epi.metadata <- amnio.epi.metadata[,c("barcode", "singleR_labels")]
write.csv(x = amnio.epi.metadata, row.names = TRUE, file = "4.1.cohort2.amnio.epi.metadata.res.0.5.20251007.csv")

####
rm(list = ls())
load("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/4.1.output/4.1.cohort2.hemato.integrated.singleR.with.inte.cluster.res.0.5.20251007.RData")
amnio.hemato.metadata <- amnio@meta.data
amnio.hemato.metadata$barcode <- rownames(amnio.hemato.metadata)
amnio.hemato.metadata <- amnio.hemato.metadata[,c("barcode", "singleR_labels")]
write.csv(x = amnio.hemato.metadata, row.names = TRUE, file = "4.1.cohort2.amnio.hemato.metadata.res.0.5.20251007.csv")
