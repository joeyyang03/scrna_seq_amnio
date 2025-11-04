rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/6.0.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 

library(dplyr)
library(readr)
library(DESeq2)
library(Seurat)
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]

celltype <- basename(input_file)
celltype <- gsub(".Robj", "", celltype)
celltype <- gsub("0.0.", "", celltype)

load(input_file)
unique(amniocytes.type@meta.data$orig.ident)

metadata <- amniocytes.type@meta.data
metadata$Trimester <- NA
metadata$Trimester[which(str_detect(metadata$orig.ident, "^GW15"))] <- "Stage.I"
metadata$Trimester[which(str_detect(metadata$orig.ident, "^GW17"))] <- "Stage.I"
metadata$Trimester[which(str_detect(metadata$orig.ident, "^GW18"))] <- "Stage.II"
metadata$Trimester[which(str_detect(metadata$orig.ident, "^GW19"))] <- "Stage.II"
metadata$Trimester[which(str_detect(metadata$orig.ident, "^GW20"))] <- "Stage.II"
metadata$Trimester[which(str_detect(metadata$orig.ident, "^GW21"))] <- "Stage.II"
metadata$Trimester[which(str_detect(metadata$orig.ident, "^GW22"))] <- "Stage.II"
metadata$Trimester[which(str_detect(metadata$orig.ident, "^GW23"))] <- "Stage.II"
metadata$Trimester[which(str_detect(metadata$orig.ident, "^GW24"))] <- "Stage.II"
metadata$Trimester[which(str_detect(metadata$orig.ident, "^GW25"))] <- "Stage.II"
metadata$Trimester[which(str_detect(metadata$orig.ident, "^GW26"))] <- "Stage.II"
metadata$Trimester[which(str_detect(metadata$orig.ident, "^GW27"))] <- "Stage.II"
metadata$Trimester[which(str_detect(metadata$orig.ident, "^GW28"))] <- "Stage.III"
metadata$Trimester[which(str_detect(metadata$orig.ident, "^GW29"))] <- "Stage.III"
metadata$Trimester[which(str_detect(metadata$orig.ident, "^GW30"))] <- "Stage.III"
amniocytes.type@meta.data <- metadata
unique(amniocytes.type@meta.data$Trimester)
table(amniocytes.type@meta.data$Trimester, amniocytes.type@meta.data$orig.ident)
amniocytes.normal.stageII <- subset(amniocytes.type,subset = Trimester %in% c("Stage.II"))
amniocytes.normal.stageII@meta.data$group <- "Normal"
###############################################################################################################################
####
load(paste0("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/5.0.output/5.0.cohort1.",celltype,".Robj"))
unique(amniocytes.type@meta.data$orig.ident)
unique(amniocytes.type@meta.data$gestation.week)
metadata <- amniocytes.type@meta.data
metadata$Trimester <- NA
metadata$Trimester[metadata$gestation.week == "GW22"] <- "Stage.II"
metadata$Trimester[metadata$gestation.week == "GW23"] <- "Stage.II"
metadata$Trimester[metadata$gestation.week == "GW24"] <- "Stage.II"
metadata$Trimester[metadata$gestation.week == "GW26"] <- "Stage.II"
metadata$Trimester[metadata$gestation.week == "GW27"] <- "Stage.II"
metadata$Trimester[metadata$gestation.week == "GW30"] <- "Stage.III"
metadata$Trimester[metadata$gestation.week == "GW31"] <- "Stage.III"
metadata$Trimester[metadata$gestation.week == "GW32"] <- "Stage.III"
amniocytes.type@meta.data <- metadata
unique(amniocytes.type@meta.data$Trimester)
table(amniocytes.type@meta.data$Trimester, amniocytes.type@meta.data$orig.ident)

amniocytes.FGR.stageII <- subset(amniocytes.type, subset = Trimester %in% c("Stage.II"))
amniocytes.FGR.stageII@meta.data$group <- "FGR"



################################################################################################################################
amniocytes.stageII <- merge(amniocytes.normal.stageII, amniocytes.FGR.stageII)

metadata <- amniocytes.stageII@meta.data
information <- metadata[,c("barcode", "group")]

merge <- amniocytes.stageII@assays$RNA@counts
merge <- as.data.frame(merge)
merge <- merge +1
dds <-DESeqDataSetFromMatrix(countData = round(merge),
                             colData = information,
                             design = ~group)
dds.1 <- DESeq(dds)
res <- results(dds.1, contrast = c("group","Normal", "FGR"))
res <- res[order(res$padj),]
res <- as.data.frame(res)
write.table(res, file = paste0("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/6.0.output/6.0.deseq2.stageII.",celltype,".20251007.txt"))


save(dds.1, file = paste0("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/6.0.output/6.0.deseq2.stageII.results.",celltype,".20251007.Robj"))

