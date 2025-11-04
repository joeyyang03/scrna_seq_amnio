rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/1.normal/7.1.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 
library(dplyr)
library("Mfuzz")
library(Seurat)
library(tibble)
library(Cairo)
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]

celltype <- basename(input_file)
celltype <- gsub(".Robj", "", celltype)
celltype <- gsub("4.0.", "", celltype)

load(input_file)
cell_metadata <- amniocytes.type@meta.data
cell_metadata$gestation.week <- cell_metadata$orig.ident
cell_metadata$gestation.week <- gsub("_.*","",cell_metadata$gestation.week)
cell_metadata <- cell_metadata[,c("barcode", "gestation.week")]

matrix <- as.data.frame(as.matrix(amniocytes.type[["RNA"]]@data))
matrix <- matrix[which(rowSums(matrix > 0.1)>=5), ]

matrix <- as.data.frame(t(matrix))
matrix$barcode <- rownames(matrix)

matrix <- merge(cell_metadata, matrix,by="barcode")
matrix <- matrix %>% select(-barcode)
matrix[1:10,1:10]

group_mean <- matrix %>%
  group_by(gestation.week) %>%
  summarise(across(where(is.numeric), ~ mean(., na.rm = TRUE), .names = "Mean_{.col}"))

group_mean <- as.data.frame(group_mean)
rownames(group_mean) <- group_mean$gestation.week
group_mean <- group_mean %>% select(-gestation.week)
group_mean[1:10,1:14]
group_mean <- as.data.frame(t(group_mean))
group_mean[1:10,1:14]
group_mean <- as.matrix(group_mean)

mfuzz_class <- new('ExpressionSet',exprs = group_mean)


mfuzz_class <- filter.NA(mfuzz_class, thres = 0.25)
mfuzz_class <- fill.NA(mfuzz_class, mode = 'mean')
mfuzz_class <- filter.std(mfuzz_class, min.std = 0)



mfuzz_class <- standardise(mfuzz_class)
set.seed(123)

cluster_num <- 15

mfuzz_cluster <- mfuzz(mfuzz_class, c = cluster_num, m = mestimate(mfuzz_class))
save.image(file = paste0("/share/home/yangjingyi/project/1.BD/up.load/1.normal/7.1.output/7.1.",celltype,".mfuzz.20251005.RData"))
