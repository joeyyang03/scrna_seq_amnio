rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/1.normal/4.0.output/")
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

load("/share/home/yangjingyi/project/1.BD/up.load/1.normal/1.0.output/1.0.amniocytes.integrated.harmony.Robj")
metadata <- amniocytes.sum1@meta.data
metadata$barcode <- rownames(metadata)
amniocytes.sum1@meta.data <- metadata

amnio.epi.metadata <- read.csv("/share/home/yangjingyi/project/1.BD/up.load/1.normal/3.2.output/3.2.amnio.epi.metadata.20251005.csv",row.names = 1)
amnio.epi.metadata <- amnio.epi.metadata[,c("barcode", "singleR_labels")]
amnio.hemato.metadata <- read.csv("/share/home/yangjingyi/project/1.BD/up.load/1.normal/3.2.output/3.2.amnio.hemato.metadata.20251005.csv",row.names = 1)
amnio.singleR <- rbind(amnio.epi.metadata, amnio.hemato.metadata)

metadata <- merge(amnio.singleR,metadata, by = "barcode")
rownames(metadata) <- metadata$barcode
amniocytes.sum1@meta.data <- metadata

cell.type <- as.data.frame(table(amnio.singleR$singleR_labels))
cell.type <- cell.type[cell.type$Freq > nrow(metadata) /100,]

amniocytes.normal <- subset(amniocytes.sum1, subset = singleR_labels %in% cell.type$Var1)

metadata.sel <- amniocytes.normal@meta.data
#####
metadata.sel$gestation.day <- NA

gestation.day <- read.csv("/share/home/yangjingyi/project/1.BD/data/rawdata.20251001/df.summary.normal.20251001.csv")
gestation.day <- gestation.day[gestation.day$Group == "Normal",]
sample <- unique(metadata.sel$orig.ident)
for (i in 1:length(sample)){
  a <- sample[[i]]
  gestation <- gestation.day[gestation.day$Sample.ID == a,]
  gestation <- gestation$Gestation.day
  metadata.sel$gestation.day[metadata.sel$orig.ident == a] <- gestation
  
}
a <- table(metadata.sel$gestation.day,metadata.sel$orig.ident)
a <- matrix(a, 
            ncol=ncol(a), 
            dimnames=dimnames(a))
a <- as.data.frame(a)

amniocytes.normal@meta.data <- metadata.sel
unique(amniocytes.normal@meta.data$gestation.day)

save(amniocytes.normal, file = "/share/home/yangjingyi/project/1.BD/up.load/1.normal/4.0.output/amniocytes.normal.Robj")

amniocytes.cell.type <- droplevels(cell.type$Var1)
for (i in 1:length(amniocytes.cell.type)){
  a <- amniocytes.cell.type[[i]]
  amniocytes.type <- subset(x = amniocytes.normal, subset = singleR_labels == a)
  save(amniocytes.type, file = paste0("/share/home/yangjingyi/project/1.BD/up.load/1.normal/4.0.output/4.0.",a,".Robj"))
  
}
