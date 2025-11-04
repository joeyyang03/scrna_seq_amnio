rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/2.model/0.0.output/")
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
load("/share/home/yangjingyi/project/1.BD/up.load/1.normal/4.0.output/amniocytes.normal.Robj")
sample.names <- unique(amniocytes.normal@meta.data$orig.ident)
strata <- gsub("_.*", "", sample.names)  
table(strata) 
df <- data.frame(sample = sample.names, strata = strata)
set.seed(123)  
selected_samples <- df %>% group_by(strata) %>% sample_n(size = min(3, n()), replace = FALSE)  
selected_samples <- selected_samples$sample

test_samples <- sample.names[!sample.names %in% selected_samples]
amniocytes.sample <- subset(amniocytes.normal, subset = orig.ident %in% selected_samples)
save(amniocytes.sample, file = paste0("/share/home/yangjingyi/project/1.BD/up.load/2.model/0.0.output/amnio.train.Robj"))

amniocytes.test <- subset(amniocytes.normal, subset = orig.ident %in% test_samples)
save(amniocytes.test, file = paste0("/share/home/yangjingyi/project/1.BD/up.load/2.model/0.0.output/amnio.test.Robj"))

amniocytes.cell.type <- unique(amniocytes.sample@meta.data$singleR_labels)
for (i in 1:length(amniocytes.cell.type)){
  a <- amniocytes.cell.type[[i]]
  amniocytes.type <- subset(x = amniocytes.sample, subset = singleR_labels == a)
  save(amniocytes.type, file = paste0("/share/home/yangjingyi/project/1.BD/up.load/2.model/0.0.output/0.0.",a,".Robj"))
  
}

