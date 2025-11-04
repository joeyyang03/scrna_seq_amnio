rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/1.normal/3.0.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 

library(Seurat)
library(ggplot2)
library(tidyverse)
library(rhdf5)
library("edgeR")
library(dplyr)
library(readr)
"%!in%" <- function(x,y)!('%in%'(x,y))
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]

cellline <- basename(input_file)
cellline <- gsub(".Robj", "", cellline)
cellline <- gsub("2.1.amnio.", "", cellline)

load(input_file)
load(paste0("/share/home/yangjingyi/project/1.BD/up.load/1.cluster/0.0.output/0.0.reference_data.",cellline,".20251004.Robj"))

seurat_list <- list(get(paste0("amnio.",cellline)),
                    get(paste0("reference_data.",cellline)))

seurat_list <- lapply(seurat_list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  dims = 1:30,
  anchor.features = 2000
)


combined_seurat <- IntegrateData(anchorset = anchors, dims = 1:30)
save.image(file = paste0("/share/home/yangjingyi/project/1.BD/up.load/1.normal/3.0.output/3.0.",cellline,".integrated.20251004.RData"))
