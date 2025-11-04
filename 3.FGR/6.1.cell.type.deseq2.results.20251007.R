rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/6.1.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 

library(dplyr)
library(readr)
library(DESeq2)
library(Seurat)
library(tidyverse)
library(ggrepel)
library(cowplot)
library(pheatmap)
library(ggplot2)
library(ggvenn)
library(reshape2)
"%!in%" <- function(x,y)!('%in%'(x,y))
folder_path <- "/share/home/yangjingyi/project/1.BD/up.load/3.FGR/6.0.output"
csv_files <- list.files(path = folder_path, pattern = "^6.0.deseq2.stage.*\\.20251007.txt$", full.names = TRUE)


matrix <- data.frame(matrix(nrow=0, ncol=6))
colnames(matrix) <- c("gene","log2FoldChange", "padj", "cell.type","sig", "stage")
for (i in 1:length(csv_files)){
  res <- read.table(csv_files[[i]],row.names = 1)
  
  stage <- gsub("6.0.deseq2.", "", basename(csv_files[[i]]))
  stage <- sub("\\..*", "", stage)
  res$stage <- stage
  
  celltype <- gsub(".20251007.txt", "", basename(csv_files[[i]]))
  celltype <- gsub("6.0.deseq2.", "", celltype)
  celltype <- sub("^[^.]*\\.", "", celltype)
  res$cell.type <- celltype
  
  
  res$sig <- NA
  res[which(res$log2FoldChange >= 1 &res$padj < 0.05), 'sig'] <- 'UP.in.Normal.group'
  res[which(res$log2FoldChange <= -1 &res$padj < 0.05), 'sig'] <- 'UP.in.FGR.group'
  res[which(abs(res$log2FoldChange) < 1 &res$padj < 0.05), 'sig'] <- 'No.significance'
  res <- res[res$sig %in% c("UP.in.Normal.group", "UP.in.FGR.group", "No.significance"),]
  
  res$gene <- rownames(res)
  res <- res[,c("gene","log2FoldChange", "padj", "cell.type", "sig", "stage")]
  matrix <- rbind(matrix, res)
}
table(matrix$cell.type)
table(matrix$stage)
able(matrix$sig,matrix$stage)
##########################################################################################################
####################################volcano.plot###########################################################
###########################################################################################################
stage <- unique(matrix$stage)  
p <- list()
for (i in 1:length(stage)){
  m <- stage[[i]]
  res <- matrix[matrix$stage == m,]
  res$p_val_adj <- res$padj
  res$p_val_adj[res$p_val_adj < 1e-300] <- 1e-300
  res <- res %>% select(-padj)
  res <- res %>% select(-sig)
  
  
  res$cell.type <- factor(res$cell.type, levels = c(
    "Stomach-Myeloid.cells",                         
    "Stomach-Erythroblasts",                         
    "Stomach-Mesothelial.cells",                     
    "Stomach-Ciliated.epithelial.cells", 
    "Kidney-Megakaryocytes", 
    "Kidney-Ureteric.bud.cells", 
    "Lung-Squamous.epithelial.cells",                
    "Lung-Bronchiolar.and.alveolar.epithelial.cells",                       
    "Intestine-Intestinal.epithelial.cells",
    "Placenta-Megakaryocytes"
  ))
  
  sig1 <- res %>%
    filter(abs(log2FoldChange) > 0.5)
  
  mycol <- c(
    "Stomach-Myeloid.cells" = "#369891", 
    "Stomach-Erythroblasts" = "#86CB66", 
    "Stomach-Mesothelial.cells" = "#91AD5A", 
    "Stomach-Ciliated.epithelial.cells" = "#007A0B", 
    "Kidney-Megakaryocytes" = "#8CCAFD", 
    "Kidney-Ureteric.bud.cells" = "#009FF6", 
    "Lung-Squamous.epithelial.cells" = "#FF662E", 
    "Lung-Bronchiolar.and.alveolar.epithelial.cells" = "#FFC763", 
    "Intestine-Intestinal.epithelial.cells" = "#FC8DCA", 
    "Placenta-Megakaryocytes" = "#704D9E"
  )
  
  p[[i]] <- ggplot() +
    geom_point(data = res,
               aes(x = -log10(p_val_adj), y = log2FoldChange),
               size = 0.8, color = 'grey') +
    geom_point(data = sig1,
               aes(x = -log10(p_val_adj), y = log2FoldChange,
                   color = cell.type),
               size = 1) +
    geom_hline(yintercept = c(-0.5, 0.5), color = "grey50", lty = 'dashed', size = 0.5) +
    geom_hline(yintercept = c(-1, 1), color = "grey50", lty = 'dashed', size = 0.5) +
    facet_grid(. ~ cell.type, scales = "free_x") +
    scale_color_manual(values = mycol) +
    ylim(min(res$log2FoldChange, na.rm = TRUE), max(res$log2FoldChange, na.rm = TRUE)) +
    theme_bw() +
    theme(
      legend.position = 'none',
      panel.grid = element_blank(),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, vjust = 0.8),
      strip.text.x = element_text(size = 10, face = 'bold')
    ) +
    labs(
      x = "-log10(adj.p)",
      y = "log2(FoldChange)",
      title = m
    )

  ggsave (p[[i]],file= paste0("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/6.1.output/6.1.deseq2.",m,".20251007.pdf"),width = 8,height = 8)
  
}

##########################################################################################################
####################################gene.number############################################################
###########################################################################################################
matrix.stageI <- matrix[matrix$stage == "stageII",]
matrix.stageI <- matrix.stageI[matrix.stageI$sig %in% c("UP.in.FGR.group", "UP.in.Normal.group"),]
matrix.stageI.plot <- as.data.frame(table(matrix.stageI$cell.type, matrix.stageI$sig))
colnames(matrix.stageI.plot) <- c("celltype", "type", "value")

matrix.stageI.plot.order <- as.data.frame(table(matrix.stageI$cell.type))
colnames(matrix.stageI.plot.order) <- c("celltype", "value")
matrix.stageI.plot.order <- matrix.stageI.plot.order[order(matrix.stageI.plot.order$value, decreasing = T),]
matrix.stageI.plot.order$celltype <- factor(matrix.stageI.plot.order$celltype, levels = matrix.stageI.plot.order$celltype)

matrix.stageI.plot$celltype <- factor(matrix.stageI.plot$celltype, levels = matrix.stageI.plot.order$celltype)

p <- ggplot(matrix.stageI.plot, aes(x = celltype, y = value, fill = type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = ifelse(value > 0, value, "")),
    position = position_stack(vjust = 0.5),  
    color = "black", size = 4          
  ) +
  scale_fill_manual(values = c("UP.in.Normal.group" = "#F8DF79", 
                               "UP.in.FGR.group" = "#B366C7")) +
  labs(x = "Cell Type",
       y = "Count") +
  theme(panel.background = element_rect(fill = 'white'), 
        axis.text.x = element_text(colour = "black", angle = 30, hjust = 1), 
        axis.text.y = element_text(colour = "black", angle = 0, hjust = 1), 
        axis.line.x = element_line(colour = "black",  size=0.5, lineend = "butt"), 
        axis.line.y = element_line(colour = "black", size=0.5))
ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/up.load/3.FGR/6.1.output/6.1.deseq2.result.stageI.Numbers.of.genes.20251007.pdf",width = 8,height = 8)


######
matrix.stageII <- matrix[matrix$stage == "stageIII",]
matrix.stageII <- matrix.stageII[matrix.stageII$sig %in% c("UP.in.FGR.group", "UP.in.Normal.group"),]
matrix.stageII.plot <- as.data.frame(table(matrix.stageII$cell.type, matrix.stageII$sig))
colnames(matrix.stageII.plot) <- c("celltype", "type", "value")

matrix.stageII.plot.order <- as.data.frame(table(matrix.stageII$cell.type))
colnames(matrix.stageII.plot.order) <- c("celltype", "value")
matrix.stageII.plot.order <- matrix.stageII.plot.order[order(matrix.stageII.plot.order$value, decreasing = T),]
matrix.stageII.plot.order$celltype <- factor(matrix.stageII.plot.order$celltype, levels = matrix.stageII.plot.order$celltype)

matrix.stageII.plot$celltype <- factor(matrix.stageII.plot$celltype, levels = matrix.stageII.plot.order$celltype)

p <- ggplot(matrix.stageII.plot, aes(x = celltype, y = value, fill = type)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = ifelse(value > 0, value, "")),
    position = position_stack(vjust = 0.5), 
    color = "black", size = 4           
  ) +
  scale_fill_manual(values = c("UP.in.Normal.group" = "#FAD797", 
                               "UP.in.FGR.group" = "#A888AF")) +
  labs(x = "Cell Type",
       y = "Count") +
  theme(panel.background = element_rect(fill = 'white'), 
        axis.text.x = element_text(colour = "black", angle = 30, hjust = 1), 
        axis.text.y = element_text(colour = "black", angle = 0, hjust = 1), 
        axis.line.x = element_line(colour = "black",  size=0.5, lineend = "butt"), 
        axis.line.y = element_line(colour = "black", size=0.5))
ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/up.load/3.FGR/6.1.output/6.1.deseq2.result.stageII.Numbers.of.genes.20251007.pdf",width = 8,height = 8)
##############################################################
######overlap##################################################
###############################################################
matrix.sel <- matrix[matrix$sig %!in% "No.significance",]
matrix.sel$sig <- sub("UP.in.Normal","Down.in.FGR",matrix.sel$sig)
unique(matrix.sel$cell.type)

celltype <- unique(matrix.sel$cell.type)
matrix <- data.frame(matrix(nrow=0, ncol=3))
colnames(matrix) <- c("gene","sig","celltype")

for (i in 1:length(celltype)){
  cell <- celltype[[i]]
  res <- matrix.sel[matrix.sel$cell.type == cell,]
  sel <- as.data.frame(table(res$gene,res$sig))
  sel <- sel[sel$Freq == 2,]
  colnames(sel) <- c("gene","sig","Fre")
  sel <- sel[,c("gene","sig")]
  if (nrow(sel) >0){
    sel$celltype <- NA
    sel$celltype <- cell
    
    matrix <- rbind(matrix,sel)
  }}
table(matrix$celltype)
table(matrix$sig)
matrix$plot <- paste0(matrix$sig,"_",matrix$celltype)
table(matrix$plot)

save(matrix, file = "/share/home/yangjingyi/project/1.BD/up.load/3.FGR/6.1.output/6.1.deseq2.results.selected.overlap.two.stages.20251007.Robj")

#############################plot
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/0.0.output/amnio.train.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/5.0.output/amniocytes.cohort1.FGR.Robj")

celltype <- unique(matrix$celltype)
matrix.summary <- data.frame(c(unique(amniocytes.sample$orig.ident),unique(amniocytes.FGR.cohort1$orig.ident)))
colnames(matrix.summary) <- "orig.ident"
for (i in 1:length(celltype)){
  a <- celltype[[i]]
  gene <- matrix[matrix$celltype == a,]$gene
  amniocytes.Normal.sel <- subset(amniocytes.sample, subset = singleR_labels == a)
  matrix.normal.sel <- as.matrix(GetAssayData(amniocytes.Normal.sel, slot = "data"))  
  matrix.normal.sel <- subset(matrix.normal.sel, rownames(matrix.normal.sel) %in% gene, drop = FALSE)
  rownames(matrix.normal.sel) <- paste0(rownames(matrix.normal.sel),"_",a)
  matrix.normal.sel <- as.data.frame(t(matrix.normal.sel))
  matrix.normal.sel$barcode <- rownames(matrix.normal.sel)
  
  
  metadata.normal.sel <- amniocytes.Normal.sel@meta.data
  metadata.normal.sel2 <- metadata.normal.sel[,c("barcode", "orig.ident")]
  
  matrix.normal.sel <- merge(matrix.normal.sel, metadata.normal.sel2,by = "barcode")
  matrix.normal.sel <- matrix.normal.sel %>% select(-barcode)
  
  matrix.normal.sel2 <- matrix.normal.sel %>%
    group_by(orig.ident) %>%
    summarise(across(everything(), mean), .groups = 'drop') 
  ###FGR
  amniocytes.FGR.sel <- subset(amniocytes.FGR.cohort1, subset = singleR_labels == a)
  matrix.FGR.sel <- as.matrix(GetAssayData(amniocytes.FGR.sel, slot = "data"))  
  matrix.FGR.sel <- subset(matrix.FGR.sel, rownames(matrix.FGR.sel) %in% gene, drop = FALSE)
  rownames(matrix.FGR.sel) <- paste0(rownames(matrix.FGR.sel),"_",a)
  matrix.FGR.sel <- as.data.frame(t(matrix.FGR.sel))
  matrix.FGR.sel$barcode <- rownames(matrix.FGR.sel)
  
  
  metadata.FGR.sel <- amniocytes.FGR.sel@meta.data
  metadata.FGR.sel2 <- metadata.FGR.sel[,c("barcode", "orig.ident")]
  
  matrix.FGR.sel <- merge(matrix.FGR.sel, metadata.FGR.sel2,by = "barcode")
  matrix.FGR.sel <- matrix.FGR.sel %>% select(-barcode)
  
  matrix.FGR.sel2 <- matrix.FGR.sel %>%
    group_by(orig.ident) %>%
    summarise(across(everything(), mean), .groups = 'drop')
  
  matrix.sel <- as.data.frame(rbind(matrix.normal.sel2,matrix.FGR.sel2))
  matrix.summary <- merge(matrix.summary,matrix.sel,by="orig.ident", all.x=T)
  print(a)
}
rownames(matrix.summary) <- matrix.summary$orig.ident
matrix.summary <- matrix.summary %>% select(-orig.ident)
matrix.summary <- matrix.summary[rownames(matrix.summary) %!in%  c("GW15_1","GW17_1"),]

annotation.1 <- amniocytes.sample@meta.data 
annotation.1 <- annotation.1[,c("orig.ident", "gestation.day")]

annotation.1 <- annotation.1 %>%
  group_by(orig.ident) %>%
  summarise(gestation.day = dplyr::first(gestation.day)) 

annotation.2 <- amniocytes.FGR.cohort1@meta.data 
annotation.2 <- annotation.2[,c("orig.ident", "gestation.day")]

annotation.2 <- annotation.2 %>%
  group_by(orig.ident) %>%
  summarise(gestation.day = dplyr::first(gestation.day)) 
annotation <- rbind(annotation.1,annotation.2)
annotation <- annotation[annotation$orig.ident %!in% c("GW15_1","GW17_1"),]
matrix.summary <- as.data.frame(t(matrix.summary))
pheatmap(matrix.summary, scale = "row", 
         color = colorRampPalette(c("navy", "white", "firebrick"))(50),
         filename = "6.1.pheatmap.overlap.genes.20251007.pdf")

