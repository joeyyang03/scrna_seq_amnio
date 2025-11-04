rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/5.0.output/")
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
######################################################################################################################
#resolution 0.5#######################################################################################################
#######################################################################################################################
load("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/1.1.output/1.1.amniocytes.integrated.harmony.cohort2.Robj")
metadata <- amniocytes.sum1@meta.data
metadata$barcode <- rownames(metadata)
amniocytes.sum1@meta.data <- metadata

amnio.epi.metadata <- read.csv("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/4.1.output/4.1.cohort2.amnio.epi.metadata.res.0.5.20251007.csv",row.names = 1)
amnio.hemato.metadata <- read.csv("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/4.1.output/4.1.cohort2.amnio.hemato.metadata.res.0.5.20251007.csv",row.names = 1)
amnio.singleR <- rbind(amnio.epi.metadata, amnio.hemato.metadata)

metadata <- merge(amnio.singleR,metadata, by = "barcode")

unique(metadata$orig.ident)
metadata$orig.ident <- gsub(".matrix", "", metadata$orig.ident)
#####
metadata$gestation.day <- NA
metadata$indicator <- NA
metadata$Frozen <- NA
metadata$Maternal.age <- NA

gestation.day <- read.csv("/share/home/yangjingyi/project/1.BD/data/rawdata.20251001/df.summary.FGR.20251008.csv")
colnames(gestation.day) <- c("Maternal.age",   "Gestation.day",  "indicator",      "Frozen",         "vali.group",     "gestation.week", "orig.ident" )
gestation.day <- gestation.day[gestation.day$vali.group == "cohort2",]

sample <- unique(metadata$orig.ident)
for (i in 1:length(sample)){
  a <- sample[[i]]
  gestation <- gestation.day[gestation.day$orig.ident == a,]
  
  metadata$gestation.day[metadata$orig.ident == a] <- gestation$Gestation.day
  metadata$indicator[metadata$orig.ident == a] <- gestation$indicator
  metadata$Frozen[metadata$orig.ident == a] <- gestation$Frozen
  metadata$Maternal.age[metadata$orig.ident == a] <- gestation$Maternal.age
  metadata$gestation.week[metadata$orig.ident == a] <- gestation$gestation.week
}
a <- table(metadata$gestation.day,metadata$orig.ident)
a <- matrix(a, 
            ncol=ncol(a), 
            dimnames=dimnames(a))
a <- as.data.frame(a)
table(is.na(metadata))

rownames(metadata) <- metadata$barcode
amniocytes.sum1@meta.data <- metadata

amniocytes.FGR.cohort2 <- amniocytes.sum1
unique(amniocytes.FGR.cohort2@meta.data$gestation.day)
save(amniocytes.FGR.cohort2, file = "/share/home/yangjingyi/project/1.BD/up.load/3.FGR/5.0.output/amniocytes.cohort2.FGR.Robj")

amniocytes.cell.type <- unique(amniocytes.FGR.cohort1@meta.data$singleR_labels)
for (i in 1:length(amniocytes.cell.type)){
  a <- amniocytes.cell.type[[i]]
  amniocytes.type <- subset(x = amniocytes.FGR.cohort1, subset = singleR_labels == a)
  save(amniocytes.type, file = paste0("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/5.0.output/5.0.cohort2.",a,".Robj"))
  
}
