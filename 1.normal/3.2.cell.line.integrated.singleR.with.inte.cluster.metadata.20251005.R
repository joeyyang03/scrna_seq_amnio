rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/1.normal/3.2.output/")
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
"%!in%" <- function(x,y)!('%in%'(x,y))

###########integrated.
load("/share/home/yangjingyi/project/1.BD/up.load/1.normal/3.1.output/3.1.epi.integrated.singleR.with.inte.cluster.20251005.RData")
metadata <- combined_seurat@meta.data
metadata$group <- NA
metadata$group[which(str_detect(metadata$orig.ident, "^GW"))] <- "amniocytes"
metadata$group[which(str_detect(metadata$orig.ident, "^Science"))] <- "Science.2020"
combined_seurat@meta.data <- metadata
p <- DimPlot(combined_seurat,group.by = "group", label = F, repel = TRUE, cols = c("#A8C58C","#467F79")) + scale_alpha_manual(values = c(0.1)) 
ggsave (p,file= "/share/home/yangjingyi/project/1.BD/up.load/1.normal/3.2.output/3.2.dimplot.amniocytes.integrated.epi.20251005.pdf",width = 8,height = 8)

amnio.epi.metadata <- amnio@meta.data
amnio.epi.metadata$barcode <- rownames(amnio.epi.metadata)
amnio.epi.metadata <- amnio.epi.metadata[,c("barcode", "singleR_labels")]
write.csv(x = amnio.epi.metadata, row.names = TRUE, file = "3.2.amnio.epi.metadata.20251005.csv")

####
amnio <- subset(combined_seurat, subset = group == "amniocytes") 
DimPlot(amnio,group.by = "group")
metadata <- amnio@meta.data
metadata$barcode <- rownames(metadata)
metadata <- merge(metadata,amnio.epi.metadata,by = "barcode")
rownames(metadata) <- metadata$barcode
amnio@meta.data <- metadata
cell.type <- as.data.frame(table(metadata$singleR_labels))
cell.type <- cell.type[cell.type$Freq > 3900,]

metadata$singleR_labels[metadata$singleR_labels %!in% cell.type$Var1] <- "Undefined"
table(metadata$singleR_labels)

amnio@meta.data <- metadata
p <- DimPlot(amnio,group.by = "singleR_labels", cols = c("#FAF3DD",  "#DCEE8C", "#A9C287","#8FC0A9", "#467F79", "#2E4F49", "Grey"))
ggsave (p,file= "/share/home/yangjingyi/project/1.BD/up.load/1.normal/3.2.output/3.2.dimplot.singleR.amniocytes.integrated.epi.20051005.pdf",width = 8,height = 8)

######################
##########hemato######
######################
rm(list = ls())
load("/share/home/yangjingyi/project/1.BD/up.load/1.normal/3.1.output/3.1.hemato.integrated.singleR.with.inte.cluster.20251005.RData")
metadata <- combined_seurat@meta.data
metadata$group <- NA
metadata$group[which(str_detect(metadata$orig.ident, "^GW"))] <- "amniocytes"
metadata$group[which(str_detect(metadata$orig.ident, "^Science"))] <- "Science.2020"
combined_seurat@meta.data <- metadata
DimPlot(combined_seurat,group.by = "group")
p <- DimPlot(combined_seurat,group.by = "group", label = F, repel = TRUE, cols = c("#CFADEA", "#7570B3")) + scale_alpha_manual(values = c(0.1))
ggsave (p,file= "/share/home/yangjingyi/project/1.BD/up.load/1.normal/3.2.output/3.2.dimplot.amniocytes.integrated.hemato.20251005.pdf",width = 8,height = 8)

amnio.hemato.metadata <- amnio@meta.data
amnio.hemato.metadata$barcode <- rownames(amnio.hemato.metadata)
amnio.hemato.metadata <- amnio.hemato.metadata[,c("barcode", "singleR_labels")]
write.csv(x = amnio.hemato.metadata, row.names = TRUE, file = "3.2.amnio.hemato.metadata.20251005.csv")


####
amnio <- subset(combined_seurat, subset = group == "amniocytes") 
DimPlot(amnio,group.by = "group")
metadata <- amnio@meta.data
metadata$barcode <- rownames(metadata)
metadata <- merge(metadata,amnio.hemato.metadata,by = "barcode")
rownames(metadata) <- metadata$barcode
amnio@meta.data <- metadata
cell.type <- as.data.frame(table(metadata$singleR_labels))
cell.type <- cell.type[cell.type$Freq > 3900,]


metadata$singleR_labels[metadata$singleR_labels %!in% cell.type$Var1] <- "Undefined"
table(metadata$singleR_labels)

amnio@meta.data <- metadata
p <- DimPlot(amnio,group.by = "singleR_labels", cols = c("#58538B", "#C25CA8","#F18791", "#F9C382", "grey"))
ggsave (p,file= "/share/home/yangjingyi/project/1.BD/up.load/1.normal/3.2.output/3.2.dimplot.singleR.amniocytes.integrated.hemato.20251005.pdf",width = 8,height = 8)
