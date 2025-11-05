rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/5.doublecheck/1.5.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 

library(Seurat)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(reshape2)
library(ggpubr)
library(SingleR)

fetal.intestine.epithelium <- readRDS("/share/home/yangjingyi/project/1.BD/data/organ.data/fetal.intestine/2020.cell/epithelium.RDS")
unique(fetal.intestine.epithelium$PCW)
unique(fetal.intestine.epithelium$Cluster)
unique(fetal.intestine.epithelium$Location)


####
#PCW have to cover to gestation week
fetal.intestine.sum <- subset(fetal.intestine.epithelium, subset = PCW %in% c(20, 17, 18, 16,19, 15, 22)) 
fetal.intestine.sum$gestation.week <- paste0("GW",fetal.intestine.sum$PCW+2)
table(fetal.intestine.sum$PCW, fetal.intestine.sum$gestation.week)

fetal.intestine.sum$group <- "Fetal.intestine"
####
load("/share/home/yangjingyi/project/1.BD/up.load/1.normal/4.0.output/amniocytes.normal.Robj")
amniocytes.intestine <- subset(amniocytes.normal, subset = singleR_labels %in% c("Intestine-Intestinal.epithelial.cells"))
metadata <- amniocytes.intestine@meta.data
metadata$group <- "amniocytes"
metadata$gestation.week <- sub("_.*", "", metadata$orig.ident)
amniocytes.intestine@meta.data <- metadata
amniocytes.intestine  <- subset(amniocytes.intestine, subset = gestation.week %in% c("GW17", "GW18", "GW19", "GW20", "GW21", "GW22", "GW24"))
table(amniocytes.intestine$gestation.week)

seurat_list <- list(amniocytes.intestine,
                    fetal.intestine.sum)


seurat_list <- lapply(seurat_list, function(x) {
  x %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
})


anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  dims = 1:30,
  anchor.features = 2000
)


combined_seurat <- IntegrateData(anchorset = anchors, dims = 1:30)
save.image(file = "/share/home/yangjingyi/project/1.BD/up.load/5.doublecheck/1.5.output/1.5.intestine.integrated.20251020.RData")
DefaultAssay(combined_seurat) <- "integrated"
combined_seurat <- ScaleData(object = combined_seurat, verbose = FALSE)
combined_seurat <- RunPCA(object = combined_seurat, npcs = 13, verbose = FALSE)
combined_seurat <- RunUMAP(object = combined_seurat, reduction = "pca", dims = 1:13)
combined_seurat <- FindNeighbors(object = combined_seurat, reduction = "pca", dims = 1:13)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)
save.image(file = "/share/home/yangjingyi/project/1.BD/up.load/5.doublecheck/1.5.output/1.5.intestine.integrated.with.inte.cluster.20251020.RData")

#load("/share/home/yangjingyi/project/1.BD/up.load/5.doublecheck/1.5.output/1.5.intestine.integrated.with.inte.cluster.20251020.RData")
p <- DimPlot(combined_seurat,group.by = "singleR_labels", label = F, repel = TRUE) + scale_alpha_manual(values = c(0.1)) 
p <- DimPlot(combined_seurat,group.by = "group", label = F, repel = TRUE) + 
  scale_alpha_manual(values = c(0.01))  +
  scale_color_manual(values = c(
    "amniocytes" = "#a5cffc",
    "Fetal.intestine" = "#f8abce"
  ))
ggsave(filename="1.5.dimplot.fetal.intestine.amnio.with.legend.20251026.pdf",p,width = 8, height = 8)

