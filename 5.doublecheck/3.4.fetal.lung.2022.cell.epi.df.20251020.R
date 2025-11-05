rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/5.doublecheck/3.4.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") #install
library(hdf5r)
library(Matrix)
library(Seurat)
library(tidyverse)
library(dplyr)
load("/share/home/yangjingyi/project/1.BD/up.load/5.doublecheck/3.4.output/3.4.lung.integrated.with.inte.cluster.20251020.RData")

matrix <- as.data.frame(as.matrix(combined_seurat[["RNA"]]@data))
####
matrix <- matrix[which(rowSums(matrix > 0.1)>=5), ]
matrix.gene <- rownames(matrix)
metadata.type <- combined_seurat@meta.data
metadata.type$barcode <- rownames(metadata.type)


######
amnio.matrix <-  data.frame(matrix(ncol = 1, nrow = 4)) 
colnames(amnio.matrix) <- c("gestation.week")
amnio.matrix$gestation.week <- c("GW17", "GW20", "GW22", "GW24")

fetal.matrix <-  data.frame(matrix(ncol = 1, nrow = 4)) 
colnames(fetal.matrix) <- c("gestation.week")
fetal.matrix$gestation.week <- c("GW17", "GW20", "GW22", "GW24")


for (m in 1:length(matrix.gene)){
  gene <- matrix.gene[[m]]
  matrix.selected <- matrix[rownames(matrix) == gene,]
  matrix.selected <- as.data.frame(t(matrix.selected))
  matrix.selected$barcode <- rownames(matrix.selected)
  matrix.merge <- merge(matrix.selected, metadata.type, by="barcode")
  
  #split the matrix
  
  ##matrix.fetal
  matrix.fetal.tmp <- matrix.merge[matrix.merge$group == "Fetal.lung",]
  matrix.fetal.tmp <- matrix.fetal.tmp[,c(gene,"gestation.week")]
  colnames(matrix.fetal.tmp) <- c("expression.level.fetal", "gestation.week")
  
  matrix.fetal.tmp <- matrix.fetal.tmp %>%
    group_by(gestation.week) %>%
    summarise(
      Mean.expr.level.fetal = mean(expression.level.fetal, na.rm = TRUE),
      gestation.week = dplyr::first(gestation.week)
    )
  if (max(matrix.fetal.tmp$Mean.expr.level.fetal) > 0.1){
    colnames(matrix.fetal.tmp) <- c("gestation.week", gene)
    fetal.matrix <- merge(fetal.matrix, matrix.fetal.tmp, by="gestation.week")
  }
  
  ####matrix.amnio
  matrix.amnio.tmp <- matrix.merge[matrix.merge$group == "amniocytes",]
  matrix.amnio.tmp <- matrix.amnio.tmp[,c(gene,"gestation.week")]
  colnames(matrix.amnio.tmp) <- c("expression.level.amnio", "gestation.week")
  matrix.amnio.tmp <- matrix.amnio.tmp[,c("expression.level.amnio", "gestation.week")]
  matrix.amnio.tmp <- matrix.amnio.tmp %>%
    group_by(gestation.week) %>%
    summarise(
      Mean.expr.level.amnio = mean(expression.level.amnio, na.rm = TRUE),
      gestation.week = dplyr::first(gestation.week)
    )
  
  if (max(matrix.amnio.tmp$Mean.expr.level.amnio) > 0.1){
    colnames(matrix.amnio.tmp) <- c("gestation.week", gene)
    amnio.matrix <- merge(amnio.matrix, matrix.amnio.tmp, by="gestation.week")
  }
  print(gene)
  
} 


write.csv(x = fetal.matrix, row.names = TRUE, file = "3.4.fetal.matrix.2020.cell.epi.20251020.csv")
write.csv(x = amnio.matrix, row.names = TRUE, file = "3.3.amnio.matrix.2020.cell.epi.20251020.csv")


