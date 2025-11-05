rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/5.doublecheck/2.3.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib")
library(hdf5r)
library(Matrix)
library(Seurat)
library(tidyverse)
library(dplyr)
load("/share/home/yangjingyi/project/1.BD/up.load/5.doublecheck/2.3.output/2.3.kidney.integrated.with.inte.cluster.20251013.RData")

matrix <- as.data.frame(as.matrix(combined_seurat[["RNA"]]@data))
####
matrix <- matrix[which(rowSums(matrix > 0.1)>=5), ]
matrix.gene <- rownames(matrix)
metadata.type <- combined_seurat@meta.data
metadata.type$barcode <- rownames(metadata.type)


######
amnio.epi.matrix <-  data.frame(matrix(ncol = 1, nrow = 2)) 
colnames(amnio.epi.matrix) <- c("gestation.week")
amnio.epi.matrix$gestation.week <- c("GW15", "GW18")

amnio.immu.matrix <-  data.frame(matrix(ncol = 1, nrow = 2)) 
colnames(amnio.immu.matrix) <- c("gestation.week")
amnio.immu.matrix$gestation.week <- c("GW15", "GW18")

fetal.epi.matrix <-  data.frame(matrix(ncol = 1, nrow = 2)) 
colnames(fetal.epi.matrix) <- c("gestation.week")
fetal.epi.matrix$gestation.week <- c("GW15", "GW18")

fetal.immu.matrix <-  data.frame(matrix(ncol = 1, nrow = 2)) 
colnames(fetal.immu.matrix) <- c("gestation.week")
fetal.immu.matrix$gestation.week <- c("GW15", "GW18")

for (m in 1:length(matrix.gene)){
  gene <- matrix.gene[[m]]
  matrix.selected <- matrix[rownames(matrix) == gene,]
  matrix.selected <- as.data.frame(t(matrix.selected))
  matrix.selected$barcode <- rownames(matrix.selected)
  matrix.merge <- merge(matrix.selected, metadata.type, by="barcode")
  
  #split the matrix
  
  ##matrix.fetal
  matrix.fetal.tmp <- matrix.merge[matrix.merge$group == "Fetal.kidney",]
  matrix.fetal.tmp <- matrix.fetal.tmp[,c(gene,"gestation.week", "compartment")]
  colnames(matrix.fetal.tmp) <- c("expression.level.fetal", "gestation.week", "compartment")
  matrix.fetal.epi.tmp <- matrix.fetal.tmp[matrix.fetal.tmp$compartment == "fetal_nephron",]
  matrix.fetal.immu.tmp <- matrix.fetal.tmp[matrix.fetal.tmp$compartment == "immune",]
  
  matrix.fetal.epi.tmp <- matrix.fetal.epi.tmp %>%
    group_by(gestation.week) %>%
    summarise(
      Mean.expr.level.fetal = mean(expression.level.fetal, na.rm = TRUE),
      gestation.week = dplyr::first(gestation.week)
    )
  if (max(matrix.fetal.epi.tmp$Mean.expr.level.fetal) > 0.1){
    colnames(matrix.fetal.epi.tmp) <- c("gestation.week", gene)
    fetal.epi.matrix <- merge(fetal.epi.matrix, matrix.fetal.epi.tmp, by="gestation.week")
  }
  
  matrix.fetal.immu.tmp <- matrix.fetal.immu.tmp %>%
    group_by(gestation.week) %>%
    summarise(
      Mean.expr.level.fetal = mean(expression.level.fetal, na.rm = TRUE),
      gestation.week = dplyr::first(gestation.week)
    )
  if (max(matrix.fetal.immu.tmp$Mean.expr.level.fetal) > 0.1){
    colnames(matrix.fetal.immu.tmp) <- c("gestation.week", gene)
    fetal.immu.matrix <- merge(fetal.immu.matrix, matrix.fetal.immu.tmp, by="gestation.week")
  }
  
  ####matrix.amnio
  matrix.amnio.tmp <- matrix.merge[matrix.merge$group == "amniocytes",]
  matrix.amnio.tmp <- matrix.amnio.tmp[,c(gene,"gestation.week", "singleR_labels")]
  colnames(matrix.amnio.tmp) <- c("expression.level.amnio", "gestation.week", "singleR_labels")
  matrix.amnio.epi.tmp <- matrix.amnio.tmp[matrix.amnio.tmp$singleR_labels == "Kidney-Ureteric.bud.cells",]
  matrix.amnio.immu.tmp <- matrix.amnio.tmp[matrix.amnio.tmp$singleR_labels == "Kidney-Megakaryocytes",]
  
  matrix.amnio.epi.tmp <- matrix.amnio.epi.tmp %>%
    group_by(gestation.week) %>%
    summarise(
      Mean.expr.level.amnio = mean(expression.level.amnio, na.rm = TRUE),
      gestation.week = dplyr::first(gestation.week)
    )
  
  if (max(matrix.amnio.epi.tmp$Mean.expr.level.amnio) > 0.1){
    colnames(matrix.amnio.epi.tmp) <- c("gestation.week", gene)
    amnio.epi.matrix <- merge(amnio.epi.matrix, matrix.amnio.epi.tmp, by="gestation.week")
  }
  
  matrix.amnio.immu.tmp <- matrix.amnio.immu.tmp %>%
    group_by(gestation.week) %>%
    summarise(
      Mean.expr.level.amnio = mean(expression.level.amnio, na.rm = TRUE),
      gestation.week = dplyr::first(gestation.week)
    )
  
  if (max(matrix.amnio.immu.tmp$Mean.expr.level.amnio) > 0.1){
    colnames(matrix.amnio.immu.tmp) <- c("gestation.week", gene)
    amnio.immu.matrix <- merge(amnio.immu.matrix, matrix.amnio.immu.tmp, by="gestation.week")
  }
  print(gene)
  
} 

write.csv(x = fetal.immu.matrix, row.names = TRUE, file = "2.3.fetal.immu.matrix.2020.cell.20251025.csv")
write.csv(x = fetal.epi.matrix, row.names = TRUE, file = "2.3.fetal.epi.matrix.2020.cell.20251025.csv")
write.csv(x = amnio.immu.matrix, row.names = TRUE, file = "2.3.amnio.immu.matrix.2020.cell.20251025.csv")
write.csv(x = amnio.epi.matrix, row.names = TRUE, file = "2.3.amnio.epi.matrix.2020.cell.20251025.csv")


