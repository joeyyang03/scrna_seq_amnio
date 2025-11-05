rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/4.model/1.0.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 


library(pROC)
library(dplyr)
library(readr)
library(Seurat)
library(tidyverse)
library(ggrepel)
library(ggplot2)
library(reshape2)
library(caret)
library(glmnet)
load("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/6.1.output/6.1.deseq2.results.selected.overlap.two.stages.20251007.Robj")

######
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/0.0.output/amnio.train.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/5.0.output/amniocytes.cohort1.FGR.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/0.0.output/amnio.test.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/5.0.output/amniocytes.cohort2.FGR.Robj")

##########################################################################################################
#######################df.for.normal######################################################################
##########################################################################################################
celltype <- unique(matrix$celltype) 
sample <- unique(amniocytes.sample$orig.ident)
df.normal <- data.frame(matrix(nrow=length(sample), ncol=0))
df.normal$orig.ident <- sample
for (i in 1:length(celltype)){
  a <- celltype[[i]]
  gene.list <- matrix[matrix$celltype == a,]$gene
  amniocytes.Normal.sel <- subset(amniocytes.sample, subset = singleR_labels == a)
  matrix.normal.sel <- as.matrix(GetAssayData(amniocytes.Normal.sel, slot = "data"))  
  matrix.normal.sel <- subset(matrix.normal.sel, rownames(matrix.normal.sel) %in% gene.list, drop = FALSE)
  matrix.normal.sel <- as.data.frame(t(matrix.normal.sel))
  colnames(matrix.normal.sel) <- paste0(colnames(matrix.normal.sel),".",a)
  matrix.normal.sel$barcode <- rownames(matrix.normal.sel)
  
  metadata.normal.sel <- amniocytes.Normal.sel@meta.data
  metadata.normal.sel2 <- metadata.normal.sel[,c("barcode", "orig.ident")]
  
  matrix.normal.sel <- merge(matrix.normal.sel, metadata.normal.sel2,by = "barcode")
  matrix.normal.sel <- matrix.normal.sel %>% select(-barcode)
  matrix.normal.sel2 <- matrix.normal.sel %>%
    group_by(orig.ident) %>%
    summarise(across(everything(), mean))
  df.normal <- merge(df.normal, matrix.normal.sel2,by="orig.ident",all=T)
  
}
df.normal$group <- "Normal"
##########################################################################################################
#######################df.for.FGR######################################################################
##########################################################################################################
celltype <- unique(matrix$celltype) 
sample <- unique(amniocytes.FGR.cohort1$orig.ident)
df.FGR <- data.frame(matrix(nrow=length(sample), ncol=0))
df.FGR$orig.ident <- sample
for (i in 1:length(celltype)){
  a <- celltype[[i]]
  gene.list <- matrix[matrix$celltype == a,]$gene
  amniocytes.FGR.sel <- subset(amniocytes.FGR.cohort1, subset = singleR_labels == a)
  matrix.FGR.sel <- as.matrix(GetAssayData(amniocytes.FGR.sel, slot = "data"))  
  matrix.FGR.sel <- subset(matrix.FGR.sel, rownames(matrix.FGR.sel) %in% gene.list, drop = FALSE)
  matrix.FGR.sel <- as.data.frame(t(matrix.FGR.sel))
  colnames(matrix.FGR.sel) <- paste0(colnames(matrix.FGR.sel),".",a)
  matrix.FGR.sel$barcode <- rownames(matrix.FGR.sel)
  
  metadata.FGR.sel <- amniocytes.FGR.sel@meta.data
  metadata.FGR.sel2 <- metadata.FGR.sel[,c("barcode", "orig.ident")]
  
  matrix.FGR.sel <- merge(matrix.FGR.sel, metadata.FGR.sel2,by = "barcode")
  matrix.FGR.sel <- matrix.FGR.sel %>% select(-barcode)
  matrix.FGR.sel2 <- matrix.FGR.sel %>%
    group_by(orig.ident) %>%
    summarise(across(everything(), mean))
  df.FGR <- merge(df.FGR, matrix.FGR.sel2,by="orig.ident",all=T)
}
df.FGR$group <- "FGR"
df <- rbind(df.normal,df.FGR)
df[is.na(df)] <- 0
#df <- df %>% select(-orig.ident)

save(df, file = "/share/home/yangjingyi/project/1.BD/up.load/4.model/1.0.output/1.0.amniocytes.normal.fgr.49.variant.df.20251007.Robj")

##########################################################################################################
#######################df.for.FGR.cohort2 and normal.test######################################################################
##########################################################################################################
celltype <- unique(matrix$celltype) 
sample <- unique(amniocytes.test$orig.ident)
df.normal.test <- data.frame(matrix(nrow=length(sample), ncol=0))
df.normal.test$orig.ident <- sample
for (i in 1:length(celltype)){
  a <- celltype[[i]]
  gene.list <- matrix[matrix$celltype == a,]$gene
  amniocytes.test.sel <- subset(amniocytes.test, subset = singleR_labels == a)
  matrix.test.sel <- as.matrix(GetAssayData(amniocytes.test.sel, slot = "data"))  
  matrix.test.sel <- subset(matrix.test.sel, rownames(matrix.test.sel) %in% gene.list, drop = FALSE)
  matrix.test.sel <- as.data.frame(t(matrix.test.sel))
  colnames(matrix.test.sel) <- paste0(colnames(matrix.test.sel),".",a)
  matrix.test.sel$barcode <- rownames(matrix.test.sel)
  
  metadata.test.sel <- amniocytes.test.sel@meta.data
  metadata.test.sel2 <- metadata.test.sel[,c("barcode", "orig.ident")]
  
  matrix.test.sel <- merge(matrix.test.sel, metadata.test.sel2,by = "barcode")
  matrix.test.sel <- matrix.test.sel %>% select(-barcode)
  matrix.test.sel2 <- matrix.test.sel %>%
    group_by(orig.ident) %>%
    summarise(across(everything(), mean))
  df.normal.test <- merge(df.normal.test, matrix.test.sel2,by="orig.ident",all=T)
}

df.normal.test[is.na(df.normal.test)] <- 0
df.normal.test$group <- "Normal"

######FGR.cohort2
celltype <- unique(matrix$celltype) 
sample <- unique(amniocytes.FGR.cohort2$orig.ident)
df.FGR.cohort2 <- data.frame(matrix(nrow=length(sample), ncol=0))
df.FGR.cohort2$orig.ident <- sample
for (i in 1:length(celltype)){
  a <- celltype[[i]]
  gene.list <- matrix[matrix$celltype == a,]$gene
  amniocytes.FGR.sel <- subset(amniocytes.FGR.cohort2, subset = singleR_labels == a)
  matrix.FGR.sel <- as.matrix(GetAssayData(amniocytes.FGR.sel, slot = "data"))  
  matrix.FGR.sel <- subset(matrix.FGR.sel, rownames(matrix.FGR.sel) %in% gene.list, drop = FALSE)
  matrix.FGR.sel <- as.data.frame(t(matrix.FGR.sel))
  colnames(matrix.FGR.sel) <- paste0(colnames(matrix.FGR.sel),".",a)
  matrix.FGR.sel$barcode <- rownames(matrix.FGR.sel)
  
  metadata.FGR.sel <- amniocytes.FGR.sel@meta.data
  metadata.FGR.sel2 <- metadata.FGR.sel[,c("barcode", "orig.ident")]
  
  matrix.FGR.sel <- merge(matrix.FGR.sel, metadata.FGR.sel2,by = "barcode")
  matrix.FGR.sel <- matrix.FGR.sel %>% select(-barcode)
  matrix.FGR.sel2 <- matrix.FGR.sel %>%
    group_by(orig.ident) %>%
    summarise(across(everything(), mean))
  df.FGR.cohort2 <- merge(df.FGR.cohort2, matrix.FGR.sel2,by="orig.ident",all=T)
}

df.FGR.cohort2[is.na(df.FGR.cohort2)] <- 0
df.FGR.cohort2$group <- NA
df.FGR.cohort2$group[df.FGR.cohort2$orig.ident == "cohort2_GW18_1"] <- "FGR"
df.FGR.cohort2$group[df.FGR.cohort2$orig.ident == "cohort2_GW18_2"] <- "FGR"
df.FGR.cohort2$group[df.FGR.cohort2$orig.ident == "cohort2_GW18_3"] <- "FGR"
df.FGR.cohort2$group[df.FGR.cohort2$orig.ident == "cohort2_GW18_4"] <- "FGR"
df.FGR.cohort2$group[df.FGR.cohort2$orig.ident == "cohort2_GW19_1"] <- "FGR"
df.FGR.cohort2$group[df.FGR.cohort2$orig.ident == "cohort2_GW19_2"] <- "FGR"
df.FGR.cohort2$group[df.FGR.cohort2$orig.ident == "cohort2_GW20_1"] <- "FGR"
df.FGR.cohort2$group[df.FGR.cohort2$orig.ident == "cohort2_GW20_2"] <- "FGR"
df.FGR.cohort2$group[df.FGR.cohort2$orig.ident == "cohort2_GW22_1"] <- "Normal"
df.FGR.cohort2$group[df.FGR.cohort2$orig.ident == "cohort2_GW22_2"] <- "Normal"
df.FGR.cohort2$group[df.FGR.cohort2$orig.ident == "cohort2_GW23_1"] <- "FGR"
df.FGR.cohort2$group[df.FGR.cohort2$orig.ident == "cohort2_GW23_2"] <- "FGR"
df.FGR.cohort2$group[df.FGR.cohort2$orig.ident == "cohort2_GW23_3"] <- "Normal"
df.FGR.cohort2$group[df.FGR.cohort2$orig.ident == "cohort2_GW24_1"] <- "Normal"
df.FGR.cohort2$group[df.FGR.cohort2$orig.ident == "cohort2_GW25_1"] <- "Normal"
df.FGR.cohort2$group[df.FGR.cohort2$orig.ident == "cohort2_GW25_2"] <- "FGR"
df.FGR.cohort2$group[df.FGR.cohort2$orig.ident == "cohort2_GW27_1"] <- "FGR"
df.FGR.cohort2$group[df.FGR.cohort2$orig.ident == "cohort2_GW32_1"] <- "FGR"
df.test <- rbind(df.normal.test,df.FGR.cohort2)
save(df.test, file = "/share/home/yangjingyi/project/1.BD/up.load/4.model/1.0.output/1.0.amniocytes.49.variant.df.test.FGR.cohort2.20251009.Robj")

