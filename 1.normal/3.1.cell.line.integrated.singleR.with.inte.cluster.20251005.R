rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/1.normal/3.1.output/")
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
cellline <- gsub(".integrated.20251004.RData", "", cellline)
cellline <- gsub("3.0.", "", cellline)

load(input_file)
# Set integrated assay as default
DefaultAssay(combined_seurat) <- "integrated"
combined_seurat <- ScaleData(object = combined_seurat, verbose = FALSE)
combined_seurat <- RunPCA(object = combined_seurat, npcs = 13, verbose = FALSE)
combined_seurat <- RunUMAP(object = combined_seurat, reduction = "pca", dims = 1:13)
combined_seurat <- FindNeighbors(object = combined_seurat, reduction = "pca", dims = 1:13)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)

ref <- subset(combined_seurat, subset = orig.ident %in% c("Science.2020"))
amnio <- subset(combined_seurat, subset = orig.ident %!in% c("Science.2020"))

#####
reference_data_celltype_barcode_dot = read.csv(file = "/share/home/yangjingyi/project/1.BD/plot.final.v4/0.reference/1.0.output/reference_172_celltype_based_1000.csv", header=T, sep=",",stringsAsFactors = F,row.names = 1)
reference_data_celltype_barcode_dot$barcode = gsub("-", ".", reference_data_celltype_barcode_dot$barcode)
reference_data_celltype_barcode_dot$celltype = reference_data_celltype_barcode_dot$Organ_cell_lineage
reference_data_celltype_barcode_dot$celltype = gsub(" ", ".", reference_data_celltype_barcode_dot$celltype)

ref.metadata <- ref@meta.data
ref.metadata$barcode <- rownames(ref.metadata)
dim(ref.metadata)

reference_data_celltype_barcode_dot <- reference_data_celltype_barcode_dot[reference_data_celltype_barcode_dot$barcode %in% ref.metadata$barcode,]
reference_data_celltype_barcode_dot <- reference_data_celltype_barcode_dot[,c("barcode","celltype")]
ref.metadata <- merge(ref.metadata,reference_data_celltype_barcode_dot,by="barcode")
ref@meta.data <- ref.metadata
table(ref$celltype)


amnio$seurat_clusters <- as.factor(amnio$seurat_clusters)
levels(amnio$seurat_clusters)

results <- SingleR(
  test = amnio@assays$integrated@data, 
  ref = ref@assays$integrated@data, 
  labels = ref$celltype, 
  clusters = amnio$seurat_clusters  # 提供聚类信息
)

###
amnio$singleR_labels <- results@listData$labels[amnio$seurat_clusters]
test <- as.data.frame(table(amnio@meta.data$singleR_labels))


########
save.image(file = paste0("/share/home/yangjingyi/project/1.BD/up.load/1.normal/3.1.output/3.1.",cellline,".integrated.singleR.with.inte.cluster.20251005.RData"))
