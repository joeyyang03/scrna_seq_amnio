rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/1.0.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 
library(Seurat)
library(ggplot2)
library(tidyverse)
library(rhdf5)
library("edgeR")
library(harmony)

########################matrix for all FGR amniocytes

folder_path <- "/share/home/yangjingyi/project/1.BD/data/rawdata.20251001/"


dir <- list.dirs(folder_path, full.names = FALSE, recursive = FALSE)
cohort1_dirs <- grep("^cohort1", dir, value = TRUE)

print(cohort1_dirs)

sample.list <- list()
for (i in 1:length(cohort1_dirs)){
  counts <- Read10X(paste0(folder_path,cohort1_dirs[[i]],"/"), unique.features = TRUE)
  sample.list[[i]] <- CreateSeuratObject(counts = counts, min.cells= 3, min.features = 200, project = cohort1_dirs[[i]])
  sample.list[[i]] <- RenameCells(sample.list[[i]], add.cell.id = cohort1_dirs[[i]])
  print(cohort1_dirs[[i]])
}
amniocytes.sum <- merge(sample.list[[1]], sample.list[2:length(sample.list)])


table(amniocytes.sum@meta.data$orig.ident)

#####
HB.gene <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")

HB_m <- match(HB.gene, rownames(amniocytes.sum@assays$RNA))
HB.gene <- rownames(amniocytes.sum@assays$RNA)[HB_m]
HB.gene <- HB.gene[!is.na(HB.gene)]
amniocytes.sum[["percent.HB"]] <- PercentageFeatureSet(amniocytes.sum, features = HB.gene, assay = "RNA")

amniocytes.sum[["percent.mt"]] <-PercentageFeatureSet(object = amniocytes.sum, pattern = "^MT-")
amniocytes.sum1 <- subset(amniocytes.sum, subset = nFeature_RNA>200 & nFeature_RNA<7000 & percent.mt<25 & percent.HB <3)
amniocytes.sum1 <- NormalizeData(object = amniocytes.sum1, normalization.method = "LogNormalize",scale.factor = 10000, verbose = FALSE)

table(amniocytes.sum1@meta.data$orig.ident)

########################################### 
amniocytes.sum1 <- FindVariableFeatures(object = amniocytes.sum1, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

save.image(file = "/share/home/yangjingyi/project/1.BD/up.load/3.FGR/1.0.output/1.0.amniocytes.normal.list.cohort1.RData")

amniocytes.sum1 <- ScaleData(amniocytes.sum1) 
amniocytes.sum1 <- RunPCA(object = amniocytes.sum1, npcs = 13, verbose = FALSE)

amniocytes.sum1 <- RunHarmony(amniocytes.sum1, reduction = "pca", group.by.vars = "orig.ident", reduction.save = "harmony")
amniocytes.sum1 <- RunUMAP(object = amniocytes.sum1, reduction = "harmony", dims = 1:13, reduction.name = "umap")

amniocytes.sum1 <- FindNeighbors(object = amniocytes.sum1, reduction = "harmony", dims = 1:13)
save(amniocytes.sum1, file = "/share/home/yangjingyi/project/1.BD/up.load/3.FGR/1.0.output/1.0.amniocytes.integrated.FindNeighbors.cohort1.Robj")

amniocytes.sum1 <- FindClusters(amniocytes.sum1, resolution = 0.4)


save.image(file = "/share/home/yangjingyi/project/1.BD/up.load/3.FGR/1.0.output/1.0.amniocytes.integrated.harmony.cohort1.RData")
save(amniocytes.sum1, file = "/share/home/yangjingyi/project/1.BD/up.load/3.FGR/1.0.output/1.0.amniocytes.integrated.harmony.cohort1.Robj")
