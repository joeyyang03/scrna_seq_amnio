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
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)
library(cowplot)
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
######Go enrichment############################################
###############################################################
gestation <- c("stageII", "stageIII")
p <- list()

for (i in 1:length(gestation)){
  ges <- gestation[[i]]
  matrix.placen.ges <- matrix[matrix$cell.type == "Placenta-Megakaryocytes" & matrix$stage == ges,]
  matrix.placen.ges$log2FoldChange <- - matrix.placen.ges$log2FoldChange
  matrix.placen.ges <- matrix.placen.ges[order(matrix.placen.ges$log2FoldChange, decreasing = TRUE),]
  matrix.placen.ges$regulation <- "normal"
  matrix.placen.ges$regulation[matrix.placen.ges$log2FoldChange < -0.5 & matrix.placen.ges$padj < 0.05] <-"down.in.fgr.group"
  matrix.placen.ges$regulation[matrix.placen.ges$log2FoldChange > 0.5 & matrix.placen.ges$padj < 0.05] <-"up.in.fgr.group"
  
  table(matrix.placen.ges$regulation)
  
  genelist <- matrix.placen.ges$log2FoldChange
  names(genelist) <- matrix.placen.ges$gene
  
  head(genelist)
  tail(genelist)
  
  egmt <- GSEA(genelist, 
               TERM2GENE = geneset, 
               pvalueCutoff = 1,  
               minGSSize = 1,
               maxGSSize = 500000)
  head(egmt@result[, 1:6])
  data <- egmt@result[, c("ID", "NES", "setSize", "pvalue")]
  data <- data[order(data$NES, decreasing = TRUE), ]
  data$ID <- factor(data$ID, levels = data$ID)
  data$xlab <- 1:nrow(data) 
  
  print(summary(data$NES))
  y_min <- floor(min(data$NES))
  y_max <- ceiling(max(data$NES))
  top_pathways <- head(as.character(data$ID), 3) 
  bottom_pathways <- tail(as.character(data$ID), 3) 
  label <- c(top_pathways, bottom_pathways)
  
  data_label <- data[data$ID %in% label, ]
  data_label$col <- c("#4292c9", "#35a153", "#f26a11", "#fe9376", "#afdd8b", "#a0c9e5")
  
  fill_color <- "#8b88b7"
  p[[i]] <- ggplot(data = data, aes(x = xlab, y = NES)) +
    geom_point(
      aes(size = -log10(pvalue)),
      alpha = 1,
      shape = 21,
      stroke = 0,
      fill = fill_color
    ) +
    scale_size_continuous(
      name = "Significance\n(-log10 p-value)",
      limits = c(0, 10),
      range = c(0, 10)
    )+
    scale_alpha_continuous(range = c(0.3, 1), name = "Significance\n(-log10 p-value)") +
    labs(
      title = paste0("Placenta in ",ges," Specific Pathway Enrichment Analysis"),
      subtitle = "Based on HALLMARK gene sets",
      x = "Hallmark gene sets", 
      y = "Normalized Enrichment Score (NES)"
    ) + 
    theme_classic(base_size = 15) + 
    scale_x_continuous(breaks = seq(0, 50, by = 10), labels = seq(0, 50, by = 10))  +
    scale_y_continuous(
      breaks = seq(y_min, y_max, by = 1),
      labels = seq(y_min, y_max, by = 1),
      limits = c(y_min, y_max),
      expand = c(0.1, 0.1)  
    ) +theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
      axis.line = element_line(color = "black", size = 0.6),  
      axis.text = element_text(face = "bold"),  
      axis.title = element_text(size = 14, face = "bold"), 
      legend.title = element_text(face = "bold"), 
      legend.position = "right",  
      legend.box = "vertical",  
      panel.grid = element_blank(),  
      panel.background = element_rect(fill = "white")  
    )+geom_text_repel(data = data_label,
                      aes(x = xlab, y = NES, label = ID),
                      size =3,
                      color = data_label$col, 
                      force =20,
                      point.padding =0.5, 
                      min.segment.length = 0, 
                      hjust = 1.2, 
                      segment.color ="grey20", 
                      segment.size = 0.3, 
                      segment.alpha = 0.8, 
                      nudge_y = -0.05) 
}

p[[1]]
ggsave(p[[1]],file = "/share/home/yangjingyi/project/1.BD/up.load/3.FGR/6.1.output/6.1.deseq2.stageII.placenta.pathway.score.20251007.pdf",width = 8, height = 8)
p[[2]]
ggsave(p[[2]],file = "/share/home/yangjingyi/project/1.BD/up.load/3.FGR/6.1.output/6.1.deseq2.stageIII.placenta.pathway.score.20251007.pdf",width = 8, height = 8)
