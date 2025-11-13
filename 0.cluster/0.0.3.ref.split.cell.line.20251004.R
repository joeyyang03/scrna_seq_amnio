rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/1.cluster/0.0.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 

library(Seurat)
library(ggplot2)
library(tidyverse)
library(rhdf5)
library("edgeR")
library(dplyr)

##############matrix.all for reference.data
reference_data_celltype_barcode = read.csv(file = "/share/home/yangjingyi/project/1.BD/up.load/1.cluster/0.0.output/reference_172_celltype_based_1000.csv", header=T, sep=",",stringsAsFactors = F, row.names = 1)
reference_data_celltype_barcode$barcode <- gsub("-",".",reference_data_celltype_barcode$barcode)
dim(reference_data_celltype_barcode)
table(reference_data_celltype_barcode$organ)

reference_data_celltype_barcode$cell.line <- NA

reference_data_celltype_barcode[grepl("Corneal.and.conjunctival.epithelial.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Epithelical.cells"
reference_data_celltype_barcode[grepl("Bronchiolar.and.alveolar.epithelial.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Epithelical.cells"
reference_data_celltype_barcode[grepl("Ciliated.epithelial.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Epithelical.cells"
reference_data_celltype_barcode[grepl("Ductal.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Epithelical.cells"
reference_data_celltype_barcode[grepl("Epicardial.fat.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Epithelical.cells"
reference_data_celltype_barcode[grepl("Mesothelial.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Epithelical.cells"
reference_data_celltype_barcode[grepl("IGFBP1_DKK1.positive.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Epithelical.cells"
reference_data_celltype_barcode[grepl("PAEP_MECOM.positive.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Epithelical.cells"
reference_data_celltype_barcode[grepl("Metanephric.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Epithelical.cells"
reference_data_celltype_barcode[grepl("Goblet.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Epithelical.cells"
reference_data_celltype_barcode[grepl("Intestinal.epithelial.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Epithelical.cells"
reference_data_celltype_barcode[grepl("Squamous.epithelial.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Epithelical.cells"
reference_data_celltype_barcode[grepl("Parietal.and.chief.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Epithelical.cells"
reference_data_celltype_barcode[grepl("Ureteric.bud.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Epithelical.cells"
reference_data_celltype_barcode[grepl("PDE1C_ACSM3.positive.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Epithelical.cells"
reference_data_celltype_barcode[grepl("Neuroendocrine.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Epithelical.cells"
reference_data_celltype_barcode[grepl("Islet.endocrine.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Epithelical.cells"
reference_data_celltype_barcode[grepl("MUC13_DMBT1.positive.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Epithelical.cells"
reference_data_celltype_barcode[grepl("Acinar.cells",  reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Epithelical.cells"

reference_data_celltype_barcode[grepl("Microglia", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Hematopoietic.cells"
reference_data_celltype_barcode[grepl("Myeloid.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Hematopoietic.cells"
reference_data_celltype_barcode[grepl("Antigen.presenting.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Hematopoietic.cells"
reference_data_celltype_barcode[grepl("Megakaryocytes", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Hematopoietic.cells"
reference_data_celltype_barcode[grepl("Lymphoid.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Hematopoietic.cells"
reference_data_celltype_barcode[grepl("Hematopoietic.stem.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Hematopoietic.cells"
reference_data_celltype_barcode[grepl("Erythroblasts", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Hematopoietic.cells"
reference_data_celltype_barcode[grepl("Thymic.epithelial.cells", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Hematopoietic.cells"
reference_data_celltype_barcode[grepl("Thymocytes", reference_data_celltype_barcode$Organ_cell_lineage),]$cell.line <- "Hematopoietic.cells"
table(reference_data_celltype_barcode$cell.line)

reference_data_celltype_barcode <- reference_data_celltype_barcode[reference_data_celltype_barcode$organ %in% c("Intestine", "Kidney", "Liver", 
                                                                                                                "Lung", "Pancreas", "Placenta", 
                                                                                                                "Stomach"),]
dim(reference_data_celltype_barcode)

reference_data_epi_barcode <- reference_data_celltype_barcode[reference_data_celltype_barcode$cell.line %in% c("Epithelical.cells"),]
reference_data_hemato_barcode <- reference_data_celltype_barcode[reference_data_celltype_barcode$cell.line %in% c("Hematopoietic.cells"),]
########################################
###prepare for reference matrix#########
########################################
reference_data1 = read.csv(file = "/share/home/yangjingyi/project/1.BD/plot.final.v4/0.reference/2.0.output/reference_selected_matrix_20250326.csv", header=T, sep=",",stringsAsFactors = F, row.names = 1)
reference_data1[1:10,1:10]

rownames(reference_data1) <- reference_data1$gene
reference_data1 <- select(reference_data1,-gene)
save(reference_data1, file = "/share/home/yangjingyi/project/1.BD/up.load/1.cluster/0.0.output/0.0.reference_data.1.20251004.Robj")
#############################################epi
reference_data.epi <- reference_data1[,colnames(reference_data1) %in% reference_data_epi_barcode$barcode]
reference_data.epi <- CreateSeuratObject(counts = reference_data.epi, min.cells = 3, min.features = 200, project = "Science.2020")
reference_data.epi <- NormalizeData(object = reference_data.epi, normalization.method = "LogNormalize",scale.factor = 10000, verbose = FALSE)
reference_data.epi <- FindVariableFeatures(object = reference_data.epi, selection.method = "vst",nfeatures = 2000, verbose = FALSE)
save(reference_data.epi, file = "/share/home/yangjingyi/project/1.BD/up.load/1.cluster/0.0.output/0.0.reference_data.epi.20251004.Robj")
#############################################hemato
reference_data.hemato <- reference_data1[,colnames(reference_data1) %in% reference_data_hemato_barcode$barcode]
reference_data.hemato <- CreateSeuratObject(counts = reference_data.hemato, min.cells = 3, min.features = 200, project = "Science.2020")
reference_data.hemato <- NormalizeData(object = reference_data.hemato, normalization.method = "LogNormalize",scale.factor = 10000, verbose = FALSE)
reference_data.hemato <- FindVariableFeatures(object = reference_data.hemato, selection.method = "vst",nfeatures = 2000, verbose = FALSE)
save(reference_data.hemato, file = "/share/home/yangjingyi/project/1.BD/up.load/1.cluster/0.0.output/0.0.reference_data.hemato.20251004.Robj")
