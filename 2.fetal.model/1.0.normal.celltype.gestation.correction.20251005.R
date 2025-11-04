rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/2.model/1.0.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 

library(dplyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]

celltype <- basename(input_file)
celltype <- gsub(".Robj", "", celltype)
celltype <- gsub("0.0.", "", celltype)

load(input_file)

matrix <- as.data.frame(as.matrix(amniocytes.type[["RNA"]]@data))
matrix <- matrix[which(rowSums(matrix > 0.1)>=5), ]
matrix.gene <- rownames(matrix)
metadata.type <- amniocytes.type@meta.data

cor.matrix <-  data.frame(matrix(ncol = 7, nrow = 0)) 
colnames(cor.matrix) <- c("cell.type", "gene", "p.value", "Pearson.correlation.coefficient", "exp.max", "exp.min", "exp.mean")



for (m in 1:length(matrix.gene)){
  gene <- matrix.gene[[m]]
  matrix.selected <- matrix[rownames(matrix) == gene,]
  matrix.selected <- as.data.frame(t(matrix.selected))
  
  
  matrix.selected$barcode <- rownames(matrix.selected)
  matrix.merge <- merge(matrix.selected, metadata.type, by="barcode")
  
  matrix.merge <- matrix.merge[,c(2,4,11)]
  colnames(matrix.merge) <- c("expression.level", "orig.ident","gestation.day")
  
  group_mean <- matrix.merge %>%
    group_by(orig.ident) %>%
    summarise(
      Mean.expr.level = mean(expression.level, na.rm = TRUE),
      gestation.day = dplyr::first(gestation.day)
    )
  
  cor <- cor.test(group_mean$gestation.day, group_mean$Mean.expr.level, method = "pearson") 
  
  group_mean <- as.data.frame(group_mean)
  rownames(group_mean) <- group_mean$orig.ident
  group_mean <- group_mean %>% select(-orig.ident)
  
  exp.max <- max(group_mean$Mean.expr.level)
  exp.min <- min(group_mean$Mean.expr.level)
  exp.mean <- mean(group_mean$Mean.expr.level)
  
  cor.matrix.selected <-  data.frame(matrix(ncol = 7, nrow = 1)) 
  colnames(cor.matrix.selected) <- c("cell.type", "gene", "p.value", "Pearson.correlation.coefficient", "exp.max", "exp.min", "exp.mean")
  
  cor.matrix.selected$cell.type <- celltype
  cor.matrix.selected$gene <- gene
  cor.matrix.selected$p.value <- cor$p.value
  cor.matrix.selected$Pearson.correlation.coefficient <- cor$estimate
  cor.matrix.selected$exp.max <- exp.max
  cor.matrix.selected$exp.min <- exp.min
  cor.matrix.selected$exp.mean <- exp.mean
  
  cor.matrix <- rbind(cor.matrix,cor.matrix.selected)
  print(gene)
  
}

write.csv(x = cor.matrix, row.names = TRUE, file = paste0("1.0.amniocytes.Pearson.correlation.coefficient.selected.",celltype,".20251005.csv"))



