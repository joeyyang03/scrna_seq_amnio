rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/5.doublecheck/3.4.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib")
library(hdf5r)
library(Matrix)
library(Seurat)
library(tidyverse)
library(dplyr)

# ==========================
read_h5_to_seurat_full <- function(file_path, sparse = TRUE) {
  library(hdf5r)
  library(Matrix)
  library(Seurat)
  
  f <- H5File$new(file_path, mode = "r")
  
  X <- f[["X"]]$read()
  if (sparse) X <- Matrix(X, sparse = TRUE)
  
  cells <- f[["obs_names"]]$read()
  genes <- f[["var_names"]]$read()
  
  rownames(X) <- genes
  colnames(X) <- cells
  
  meta_names <- names(f)[grepl("^obs_", names(f)) & names(f) != "obs_names"]
  meta_list <- lapply(meta_names, function(n) f[[n]]$read())
  names(meta_list) <- gsub("^obs_", "", meta_names)  # 去掉前缀
  meta <- as.data.frame(meta_list, row.names = cells)
  
  f$close_all()
  
  seurat_obj <- CreateSeuratObject(counts = X, meta.data = meta)
  return(seurat_obj)
}

# ==========================
# ==========================

dir <- "/share/home/yangjingyi/project/1.BD/data/organ.data/fetal.lung/2022.cell/C1filtered_matrix.h5"
  
fetal.lung.sum <- read_h5_to_seurat_full(dir)
fetal.lung.sum$label <- "Epithelium(no.cilium)"
fetal.lung.sum$group <- "Fetal.lung"
fetal.lung.sum$gestation.week <- paste0("GW", as.numeric(fetal.lung.sum$stage) + 2)
table(fetal.lung.sum$gestation.week)
fetal.lung.sum <- subset(fetal.lung.sum, subset = gestation.week %in% c("GW17", "GW20", "GW22", "GW24"))
####
load("/share/home/yangjingyi/project/1.BD/up.load/1.normal/4.0.output/amniocytes.normal.Robj")
amniocytes.lung <- subset(amniocytes.normal, subset = singleR_labels %in% c("Lung-Squamous.epithelial.cells",
                                                                            "Lung-Bronchiolar.and.alveolar.epithelial.cells"))
metadata <- amniocytes.lung@meta.data
metadata$group <- "amniocytes"
metadata$gestation.week <- sub("_.*", "", metadata$orig.ident)

amniocytes.lung@meta.data <- metadata

amniocytes.lung  <- subset(amniocytes.lung, subset = gestation.week %in% c("GW17", "GW20", "GW22", "GW24"))
seurat_list <- list(amniocytes.lung,
                    fetal.lung.sum)

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
save.image(file = "/share/home/yangjingyi/project/1.BD/up.load/5.doublecheck/3.4.output/3.4.lung.integrated.20251020.RData")
DefaultAssay(combined_seurat) <- "integrated"
combined_seurat <- ScaleData(object = combined_seurat, verbose = FALSE)
combined_seurat <- RunPCA(object = combined_seurat, npcs = 13, verbose = FALSE)
combined_seurat <- RunUMAP(object = combined_seurat, reduction = "pca", dims = 1:13)
combined_seurat <- FindNeighbors(object = combined_seurat, reduction = "pca", dims = 1:13)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)
save.image(file = "/share/home/yangjingyi/project/1.BD/up.load/5.doublecheck/3.4.output/3.4.lung.integrated.with.inte.cluster.20251020.RData")


#load("/share/home/yangjingyi/project/1.BD/up.load/5.doublecheck/3.4.output/3.4.lung.integrated.with.inte.cluster.20251020.RData")
p <- DimPlot(combined_seurat,group.by = "group", label = F, repel = TRUE) + 
  scale_alpha_manual(values = c(0.1))  +
  scale_color_manual(values = c(
    "amniocytes" = "#a5cffc",
    "Fetal.lung" = "#f8abce"
  ))
ggsave(filename="3.4.dimplot.fetal.lung.amnio.with.legend.20251026.pdf",p,width = 8, height = 8)
