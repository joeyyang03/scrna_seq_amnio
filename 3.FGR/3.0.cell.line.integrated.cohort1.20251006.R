rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/3.0.output/")
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
cellline <- gsub("2.0.cohort1.amnio.", "", cellline)

load(input_file)

load("/share/home/yangjingyi/project/1.BD/up.load/1.normal/4.0.output/amniocytes.normal.Robj")
metadata <- amniocytes.normal@meta.data
metadata$cell.line <- NA
metadata[metadata$singleR_labels %in% c("Kidney-Megakaryocytes", "Placenta-Megakaryocytes", 
                                        "Stomach-Myeloid.cells", "Stomach-Erythroblasts"),]$cell.line <- "Hematopoietic.cells"
metadata[metadata$singleR_labels %in% c("Lung-Squamous.epithelial.cells",
                                        "Kidney-Ureteric.bud.cells",
                                        "Stomach-Ciliated.epithelial.cells",
                                        "Stomach-Mesothelial.cells",
                                        "Intestine-Intestinal.epithelial.cells",
                                        "Lung-Bronchiolar.and.alveolar.epithelial.cells"),]$cell.line <- "Epithelial.cells"
metadata$group <- "amniocytes"
amniocytes.normal@meta.data <- metadata
normal.epi <- subset(amniocytes.normal, subset = cell.line == "Epithelial.cells")
normal.hemato <- subset(amniocytes.normal, subset = cell.line == "Hematopoietic.cells")

seurat_list <- list(get(paste0("amnio.",cellline)),
                    get(paste0("normal.",cellline)))

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
save.image(file = paste0("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/3.0.output/3.0.cohort1.",cellline,".integrated.20251006.RData"))
