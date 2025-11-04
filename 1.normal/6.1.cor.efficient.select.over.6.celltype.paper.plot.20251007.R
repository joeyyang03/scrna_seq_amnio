rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/1.normal/6.1.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") #install

library(dplyr)
library(tidyr)
library(readr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(gg.gap)

gene.sel.up.over.6 <- read.table("/share/home/yangjingyi/project/1.BD/up.load/1.normal/6.1.output/6.1.AIGs.6.cell.types.20251007.txt")#284
gene.sel.down.over.6 <- read.table("/share/home/yangjingyi/project/1.BD/up.load/1.normal/6.1.output/6.1.ADGs.6.cell.types.20250828.txt")#43

################################################################################################
####lung########################################################################################
################################################################################################
lung.gene <- c("PSCA", "CDC20B", "AREG", "FOXJ1", "ALOX15", "IL33", "IGFBP4", "IGFBP6", "KISS1", "SCGB3A2", "NXPH4",
               "SPDEF", "SCGB3A1", "MUC16", "SCGB1A1", "LTF", "KRT4", "APOA1", "SOX9", "TESC", "SFTPA1", "SFTPC", "AGER",
               "CFHR1", "CYTL1", "UCHL1", "ASCL1", "MUC5AC", "GHRL")#cell A human fetal lung cell atlas uncovers proximal-distal gradients of differentiation and key regulators of epithelial fates
table(lung.gene %in% gene.sel.up.over.6$V1)
  
gene.sel.up.over.6[gene.sel.up.over.6$V1 %in% lung.gene,]

table(lung.gene %in% gene.sel.down.over.6$V1)



################################################################################################
####intestine########################################################################################
################################################################################################
intestine.gene <- read.csv("/share/home/yangjingyi/project/1.BD/plot.final.v7/1.normal.cell.line.filt.organ/11.2.output/fetal.intestine.20250831.csv", header = T)
table(intestine.gene$genes %in% gene.sel.up.over.6$V1)

intestine.gene[intestine.gene$genes %in% gene.sel.up.over.6$V1,]


table(intestine.gene$genes %in% gene.sel.down.over.6$V1)

intestine.gene[intestine.gene$genes %in% gene.sel.down.over.6$V1,]


################################################################################################
####stomach########################################################################################
################################################################################################
stomach.gene <- read.csv("/share/home/yangjingyi/project/1.BD/plot.final.v7/1.normal.cell.line.filt.organ/11.2.output/11.2.fetal.stomach.20250831.csv", header = T)
table(stomach.gene$gene %in% gene.sel.up.over.6$V1)

stomach.gene[stomach.gene$gene %in% gene.sel.up.over.6$V1,]

table(stomach.gene$gene %in% gene.sel.down.over.6$V1)

stomach.gene[stomach.gene$gene %in% gene.sel.down.over.6$V1,]
