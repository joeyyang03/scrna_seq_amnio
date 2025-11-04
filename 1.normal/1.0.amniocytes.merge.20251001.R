rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/1.normal/1.0.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 

library(Seurat)
library(ggplot2)
library(tidyverse)
library(rhdf5)
library("edgeR")
library(harmony)

#######################matrix for all amniocytes
samples <- read.csv("/share/home/yangjingyi/project/1.BD/data/rawdata.20251001/df.summary.normal.20251001.csv")
samples <- samples[samples$Group == "Normal",]
samples <- samples$Sample.ID
sample.list <- list()
dir <- paste0("/share/home/yangjingyi/project/1.BD/data/rawdata.20251001/", samples,".matrix/")
for (i in 1:length(dir)){
  counts <- Read10X(dir[[i]], unique.features = TRUE)
  sample.list[[i]] <- CreateSeuratObject(counts = counts, min.cells= 3, min.features = 200, project = samples[[i]])
  sample.list[[i]] <- RenameCells(sample.list[[i]], add.cell.id = samples[[i]])
  
}
amniocytes.sum <- merge(sample.list[[1]], sample.list[2:length(sample.list)])

## Print the overview of cell counts per sample
table(amniocytes.sum@meta.data$orig.ident)
###################
amniocytes.sum[["percent.mt"]] <-PercentageFeatureSet(object = amniocytes.sum, pattern = "^MT-")

HB.gene <- c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
HB_m <- match(HB.gene, rownames(amniocytes.sum@assays$RNA))
HB.gene <- rownames(amniocytes.sum@assays$RNA)[HB_m]
HB.gene <- HB.gene[!is.na(HB.gene)]
amniocytes.sum[["percent.HB"]] <- PercentageFeatureSet(amniocytes.sum, features = HB.gene, assay = "RNA")

amniocytes.sum1 <- subset(amniocytes.sum, subset = nFeature_RNA>200 & nFeature_RNA<7000 & percent.mt<25 & percent.HB <3)
amniocytes.sum1 <- NormalizeData(object = amniocytes.sum1, normalization.method = "LogNormalize",scale.factor = 10000, verbose = FALSE)

table(amniocytes.sum1@meta.data$orig.ident)
########################################### yingshe
amniocytes.sum1 <- FindVariableFeatures(object = amniocytes.sum1, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

save.image(file = "/share/home/yangjingyi/project/1.BD/up.load/1.normal/1.0.output/1.0.amniocytes.normal.list.RData")

amniocytes.sum1 <- ScaleData(amniocytes.sum1) 
amniocytes.sum1 <- RunPCA(object = amniocytes.sum1, npcs = 13, verbose = FALSE)

amniocytes.sum1 <- RunHarmony(amniocytes.sum1, reduction = "pca", group.by.vars = "orig.ident", reduction.save = "harmony")
amniocytes.sum1 <- RunUMAP(object = amniocytes.sum1, reduction = "harmony", dims = 1:13, reduction.name = "umap")

amniocytes.sum1 <- FindNeighbors(object = amniocytes.sum1, reduction = "harmony", dims = 1:13)
save(amniocytes.sum1, file = "/share/home/yangjingyi/project/1.BD/up.load/1.normal/1.0.output/1.0.amniocytes.integrated.FindNeighbors.Robj")

amniocytes.sum1 <- FindClusters(amniocytes.sum1, resolution = 0.4)


save.image(file = "/share/home/yangjingyi/project/1.BD/up.load/1.normal/1.0.output/1.0.amniocytes.integrated.harmony.RData")
save(amniocytes.sum1, file = "/share/home/yangjingyi/project/1.BD/up.load/1.normal/1.0.output/1.0.amniocytes.integrated.harmony.Robj")
