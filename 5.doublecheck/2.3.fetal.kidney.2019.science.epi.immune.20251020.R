rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/5.doublecheck/2.3.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib")
library(hdf5r)
library(Matrix)
library(Seurat)
library(tidyverse)
library(SingleR)


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
  names(meta_list) <- gsub("^obs_", "", meta_names)  
  meta <- as.data.frame(meta_list, row.names = cells)
  
  f$close_all()
  

  seurat_obj <- CreateSeuratObject(counts = X, meta.data = meta)
  return(seurat_obj)
}

# ==========================
file_path <- "/share/home/yangjingyi/project/1.BD/data/organ.data/fetal.kidney/2019.science/Fetal_full_v3.matrix.h5"
f <- H5File$new(file_path, mode = "r")
f$ls(recursive = TRUE)

fetal.kidney <- read_h5_to_seurat_full(file_path)
unique(fetal.kidney$Age)
unique(fetal.kidney$celltype)
unique(fetal.kidney$Short_Sample)
unique(fetal.kidney$compartment)
fetal.kidney <- subset(fetal.kidney, subset = Age %in% c("16+0", "13+6"))
fetal.kidney$gestation.week <- paste0("GW",as.numeric(sub("\\+.*", "", fetal.kidney$Age))+2)
table(fetal.kidney$gestation.week,fetal.kidney$Age)
fetal.kidney$group <- "Fetal.kidney"
fetal.kidney <- subset(fetal.kidney, subset = compartment %in% c("fetal_nephron", "immune"))
####
load("/share/home/yangjingyi/project/1.BD/up.load/1.normal/4.0.output/amniocytes.normal.Robj")
amniocytes.kidney <- subset(amniocytes.normal, subset = singleR_labels %in% c("Kidney-Megakaryocytes", 
                                                                              "Kidney-Ureteric.bud.cells"))
metadata <- amniocytes.kidney@meta.data
metadata$group <- "amniocytes"
metadata$gestation.week <- sub("_.*", "", metadata$orig.ident)

amniocytes.kidney@meta.data <- metadata

amniocytes.kidney  <- subset(amniocytes.kidney, subset = gestation.week %in% c("GW15", "GW18"))


seurat_list <- list(amniocytes.kidney,
                    fetal.kidney)

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
save.image(file = "/share/home/yangjingyi/project/1.BD/up.load/5.doublecheck/2.3.output/2.3.kidney.integrated.20251013.RData")
########
DefaultAssay(combined_seurat) <- "integrated"
combined_seurat <- ScaleData(object = combined_seurat, verbose = FALSE)
combined_seurat <- RunPCA(object = combined_seurat, npcs = 13, verbose = FALSE)
combined_seurat <- RunUMAP(object = combined_seurat, reduction = "pca", dims = 1:13)
combined_seurat <- FindNeighbors(object = combined_seurat, reduction = "pca", dims = 1:13)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)
save.image(file = "/share/home/yangjingyi/project/1.BD/up.load/5.doublecheck/2.3.output/2.3.kidney.integrated.with.inte.cluster.20251013.RData")


#load("/share/home/yangjingyi/project/1.BD/up.load/5.doublecheck/2.3.output/2.3.kidney.integrated.with.inte.cluster.20251013.RData")
p <- DimPlot(combined_seurat,group.by = "group", label = F, repel = TRUE) + 
  scale_alpha_manual(values = c(0.1))  +
  scale_color_manual(values = c(
    "amniocytes" = "#a5cffc",
    "Fetal.kidney" = "#f8abce"
  ))
ggsave(filename="2.3.dimplot.fetal.kidney.amnio.with.legend.20251026.pdf",p,width = 8, height = 8)


