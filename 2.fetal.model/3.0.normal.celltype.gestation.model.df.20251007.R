rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/2.model/3.0.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib")

library(dplyr)
library(readr)
library(Seurat)
library(ggplot2)
library(pheatmap)
library(cowplot)
####selection
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/1.0.output/1.0.amniocytes.with.celltype.cor.exp.matrix.20251007.Robj")
dim(cor.matrix)

cor.matrix <- cor.matrix[cor.matrix$exp.mean > 0.1, ]
table(cor.matrix$cell.type)

cor.matrix.sel <- cor.matrix[cor.matrix$p.value < 0.05, ]
table(cor.matrix.sel$cell.type)
cor.matrix.sel$organ <- cor.matrix.sel$cell.type
cor.matrix.sel$organ <- gsub("-.*", "", cor.matrix.sel$organ)
table(cor.matrix.sel$organ)

cor.matrix.sel2 <- cor.matrix.sel %>%
  group_by(organ)%>%
  top_n(n = 30, 
        wt = abs(Pearson.correlation.coefficient))
table(cor.matrix.sel2$cell.type)


cor.gene <- unique(cor.matrix.sel2$gene)

cor.cell.type <- unique(cor.matrix.sel2$cell.type)

load("/share/home/yangjingyi/project/1.BD/up.load/2.model/0.0.output/amnio.train.Robj")
matrix <- as.data.frame(as.matrix(amniocytes.sample[["RNA"]]@data))#239591 variables
matrix.sel <- matrix[rownames(matrix) %in% cor.gene,]
matrix.sel[1:10,1:10]

metadata <- amniocytes.sample@meta.data
unique(metadata$orig.ident)
metadata <- metadata[,c("barcode", "orig.ident", "singleR_labels")]
colnames(metadata) <- c("barcode", "orig.ident", "cell.type")
metadata <- metadata[metadata$cell.type %in% cor.cell.type,]


df <- data.frame(matrix(nrow=length(unique(metadata$orig.ident)),ncol=0))
df$orig.ident <- unique(amniocytes.sample@meta.data$orig.ident)

for (i in 1:length(cor.cell.type)){
  celltype.sel <- cor.cell.type[[i]]
  gene.sel <- cor.matrix.sel2[cor.matrix.sel2$cell.type == celltype.sel,]
  gene.sel <- gene.sel$gene
  matrix.sel <- matrix[rownames(matrix) %in% gene.sel,]
  
  metadata.sel <- metadata[metadata$cell.type == celltype.sel,]
  matrix.sel <- matrix.sel[,colnames(matrix.sel) %in% metadata.sel$barcode]
  
  
  for (m in 1:length(gene.sel)){
    gene <- gene.sel[[m]]
    matrix.selected <- matrix.sel[rownames(matrix.sel) == gene,]
    matrix.selected <- as.data.frame(t(matrix.selected))
    colnames(matrix.selected) <- c("Expression.level")

    
    matrix.selected$barcode <- rownames(matrix.selected)
    matrix.merge <- merge(matrix.selected, metadata.sel, by="barcode")

    group_mean <- matrix.merge %>%
      group_by(orig.ident) %>%
      summarise_at(vars(Expression.level),
                   list(Mean.expr.level = mean))
    
    colnames(group_mean) <- c("orig.ident",paste0(gene,"_",celltype.sel))
    
    df <- merge(df,group_mean,by="orig.ident", all = T)
    
    print(gene)
  }
}
dim(df)
save(df, file = "/share/home/yangjingyi/project/1.BD/up.load/2.model/3.0.output/3.0.amniocytes.with.celltype.cor.top30.exp.0.1.df.2025007.Robj")

#############################################################amnio.test.df###################
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/0.0.output/amnio.test.Robj")
amniocytes.cell.type <- unique(amniocytes.test@meta.data$singleR_labels)

matrix <- as.data.frame(as.matrix(amniocytes.test[["RNA"]]@data))
matrix <- matrix[rownames(matrix) %in% cor.gene,]

metadata <- amniocytes.test@meta.data
unique(metadata$orig.ident)
metadata <- metadata[,c("barcode", "orig.ident", "singleR_labels")]
colnames(metadata) <- c("barcode", "orig.ident", "cell.type")
metadata <- metadata[metadata$cell.type %in% cor.cell.type,]

sample <- unique(metadata$orig.ident)
df.test <- data.frame(matrix(nrow=length(sample),ncol=0))
df.test$orig.ident <- unique(amniocytes.test@meta.data$orig.ident)


for (i in 1:length(cor.cell.type)){
  celltype.sel <- cor.cell.type[[i]]
  gene.sel <- cor.matrix.sel2[cor.matrix.sel2$cell.type == celltype.sel,]
  gene.sel <- gene.sel$gene
  matrix.sel <- matrix[rownames(matrix) %in% gene.sel,]
  
  metadata.sel <- metadata[metadata$cell.type == celltype.sel,]
  matrix.sel <- matrix.sel[,colnames(matrix.sel) %in% metadata.sel$barcode]
  
  
  for (m in 1:length(gene.sel)){
    gene <- gene.sel[[m]]
    matrix.selected <- matrix.sel[rownames(matrix.sel) == gene,]
    matrix.selected <- as.data.frame(t(matrix.selected))
    colnames(matrix.selected) <- c("Expression.level")
    
    matrix.selected$barcode <- rownames(matrix.selected)
    matrix.merge <- merge(matrix.selected, metadata.sel, by="barcode")
    
    group_mean <- matrix.merge %>%
      group_by(orig.ident) %>%
      summarise_at(vars(Expression.level),
                   list(Mean.expr.level = mean))
    
    colnames(group_mean) <- c("orig.ident",paste0(gene,"_",celltype.sel))
    
    df.test <- merge(df.test,group_mean,by="orig.ident", all = T)
    
    print(gene)
  }
}
save(df.test, file = "/share/home/yangjingyi/project/1.BD/up.load/2.model/3.0.output/3.0.amniocytes.test.with.celltype.cor.top30.exp.0.1.df.20251007.Robj")
#############################################################amnio.FGR.cohort1###################
load("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/5.0.output/amniocytes.cohort1.FGR.Robj")
amniocytes.cell.type <- unique(amniocytes.FGR.cohort1@meta.data$singleR_labels)
 

matrix <- as.data.frame(as.matrix(amniocytes.FGR.cohort1[["RNA"]]@data))
matrix <- matrix[rownames(matrix) %in% cor.gene,]

metadata <- amniocytes.FGR.cohort1@meta.data
unique(metadata$orig.ident)
metadata <- metadata[,c("barcode", "orig.ident", "singleR_labels")]
colnames(metadata) <- c("barcode", "orig.ident", "cell.type")
metadata <- metadata[metadata$cell.type %in% cor.cell.type,]

sample <- unique(metadata$orig.ident)
df.FGR.cohort1 <- data.frame(matrix(nrow=length(sample),ncol=0))
df.FGR.cohort1$orig.ident <- unique(amniocytes.FGR.cohort1@meta.data$orig.ident)


for (i in 1:length(cor.cell.type)){
  celltype.sel <- cor.cell.type[[i]]
  gene.sel <- cor.matrix.sel2[cor.matrix.sel2$cell.type == celltype.sel,]
  gene.sel <- gene.sel$gene
  matrix.sel <- matrix[rownames(matrix) %in% gene.sel,]
  
  metadata.sel <- metadata[metadata$cell.type == celltype.sel,]
  matrix.sel <- matrix.sel[,colnames(matrix.sel) %in% metadata.sel$barcode]
  
  
  for (m in 1:length(gene.sel)){
    gene <- gene.sel[[m]]
    matrix.selected <- matrix.sel[rownames(matrix.sel) == gene,]
    matrix.selected <- as.data.frame(t(matrix.selected))
    colnames(matrix.selected) <- c("Expression.level")
    
    matrix.selected$barcode <- rownames(matrix.selected)
    matrix.merge <- merge(matrix.selected, metadata.sel, by="barcode")
    
    group_mean <- matrix.merge %>%
      group_by(orig.ident) %>%
      summarise_at(vars(Expression.level),
                   list(Mean.expr.level = mean))
    
    colnames(group_mean) <- c("orig.ident",paste0(gene,"_",celltype.sel))
    
    df.FGR.cohort1 <- merge(df.FGR.cohort1,group_mean,by="orig.ident", all = T)
    
    print(gene)
  }
}
save(df.FGR.cohort1, file = "/share/home/yangjingyi/project/1.BD/up.load/2.model/3.0.output/3.0.amniocytes.FGR.cohort1.with.celltype.cor.top30.exp.0.1.df.20251007.Robj")
