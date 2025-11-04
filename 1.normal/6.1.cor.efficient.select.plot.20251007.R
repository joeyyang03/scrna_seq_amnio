rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/1.normal/6.1.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib")

library(dplyr)
library(tidyr)
library(readr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(gg.gap)
####selection
folder_path <- "/share/home/yangjingyi/project/1.BD/up.load/1.normal/6.0.output"
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

cor.matrix <-  data.frame(matrix(ncol = 7, nrow = 0)) 
colnames(cor.matrix) <- c("cell.type", "gene", "p.value", "Pearson.correlation.coefficient", "exp.max","exp.min", "exp.mean")
for (i in 1:length(csv_files)){
  res <- read.csv(csv_files[[i]],row.names = 1)
  res <- res[res$p.value < 0.05,]
  cor.matrix <- rbind(cor.matrix,res)
  
}
dim(cor.matrix)

cor.matrix <- cor.matrix[cor.matrix$exp.mean > 0.1, ]
table(cor.matrix$cell.type)


cor.matrix.sel <- cor.matrix[cor.matrix$p.value < 0.05, ]
table(cor.matrix.sel$cell.type)

####################################################################################################################################################
####################################for plot########################################################################################################
####################################################################################################################################################
cor.matrix.plot <- cor.matrix.sel
cor.matrix.plot <- cor.matrix.plot[abs(cor.matrix.plot$Pearson.correlation.coefficient) > 0.5,]

cor.matrix.plot$Correlation[cor.matrix.plot$Pearson.correlation.coefficient > 0.5] <- "AIG"
cor.matrix.plot$Correlation[cor.matrix.plot$Pearson.correlation.coefficient < -0.5] <- "ADG"
table(cor.matrix.plot$Correlation)

cor.plot <- as.data.frame(table(cor.matrix.plot$cell.type, cor.matrix.plot$Correlation))
colnames(cor.plot) <- c("celltype", "Correlation", "value")

cor.plot.order <- as.data.frame(table(cor.matrix.plot$cell.type))
colnames(cor.plot.order) <- c("celltype", "value")
cor.plot.order <- cor.plot.order[order(cor.plot.order$value, decreasing = T),]
cor.plot.order$celltype <- factor(cor.plot.order$celltype, levels = cor.plot.order$celltype)

cor.plot$celltype <- factor(cor.plot$celltype, levels = cor.plot.order$celltype)

p <- ggplot(cor.plot, aes(x = celltype, y = value, fill = Correlation)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("AIG" = "#F58A8E", 
                               "ADG" = "#8A90D0"))+
  labs(x = "Cell Type",
       y = "Count") +
  theme(panel.background = element_rect(fill = 'white'), 
        axis.text.x = element_text(colour = "black",angle = 30, hjust = 1), 
        axis.text.y = element_text(colour = "black",angle = 0, hjust = 1), 
        axis.line.x = element_line(colour = "black",  size=0.5, lineend = "butt"), 
        axis.line.y = element_line(colour = "black", size=0.5)) +
  geom_text(
    mapping = aes(x = reorder(celltype, -value), y = value, label = value),
    stat = "identity",
    vjust = -0.8, hjust = 0.5)
ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/up.load/1.normal/6.1.output/6.1.cor.0.5.exp.mean.0.1.Numbers.of.genes.without.cut.20251007.pdf",width = 8,height = 8)


table(abs(cor.matrix.plot$Pearson.correlation.coefficient) > 0.5)

table(cor.matrix.plot$p.value < 0.05)

unique(cor.matrix.plot$cell.type)



######################select####################################
cor.matrix.up <- cor.matrix.sel[cor.matrix.sel$Pearson.correlation.coefficient > 0.5,]
cor.matrix.down <- cor.matrix.sel[cor.matrix.sel$Pearson.correlation.coefficient < -0.5,]

################################################################
#####up for plot################################################
################################################################
cor.matrix.up.data <- as.data.frame(table(cor.matrix.up$gene))
colnames(cor.matrix.up.data) <- c("Gene", "freq")


for (i in 1:5){
  gene.sel <- cor.matrix.up.data[cor.matrix.up.data$freq == i,]
  gene.num <- nrow(gene.sel)
  cell_order <- as.character(cor.plot.order$celltype)
  
  
  matrix.plot.fre <- cor.matrix.up[cor.matrix.up$gene %in% gene.sel$Gene,]
  matrix.plot.fre <- matrix.plot.fre[,c("gene",  "cell.type", "Pearson.correlation.coefficient")]
  matrix.plot.fre <- matrix.plot.fre %>%
    mutate(cell.type = factor(cell.type, levels = cell_order)) %>%
    arrange(cell.type, gene)
  
  length(unique(matrix.plot.fre$cell.type))
  
  matrix.plot.fre <- matrix.plot.fre %>%
    pivot_wider(names_from = cell.type,
                values_from = Pearson.correlation.coefficient) %>%
    as.data.frame()
  rownames(matrix.plot.fre) <- matrix.plot.fre$gene
  matrix.plot.fre <- matrix.plot.fre[ , -1] 
  
  missing_cols <- setdiff(cell_order, colnames(matrix.plot.fre))
  if(length(missing_cols) > 0){
    matrix.plot.fre[ , missing_cols] <- 0
  }
  
  matrix.plot.fre <- matrix.plot.fre[ , cell_order, drop = FALSE]
  matrix.plot.fre[is.na(matrix.plot.fre)] <- 0
  
  breaksList <- c(seq(0, 0.8, length.out=30)[-1],
                  seq(0.8, 1, length.out=10)[-1])
  my_palette <- colorRampPalette(c("white","#FCE8DD","#E9212C"))(length(breaksList)-1)
  pheatmap(matrix.plot.fre,
           cluster_rows = F,
           cluster_cols = F,
           color = my_palette,
           border_color = NA,
           breaks = breaksList,
           main = paste0(i," cell types, ", gene.num," AIGs"), 
           filename = paste0("6.1.pheatmap.genes.ovelap.",i," cell types, ", gene.num," AIGs.pdf"))
  
}

gene.sel <- cor.matrix.up.data[cor.matrix.up.data$freq > 5,]
write.table(gene.sel$Gene, file = "/share/home/yangjingyi/project/1.BD/up.load/1.normal/6.1.output/6.1.AIGs.6.cell.types.20251007.txt", quote = F,row.names = F,col.names = F)

gene.num <- nrow(gene.sel)

cell_order <- as.character(cor.plot.order$celltype)


matrix.plot.fre <- cor.matrix.up[cor.matrix.up$gene %in% gene.sel$Gene,]
matrix.plot.fre <- matrix.plot.fre[,c("gene",  "cell.type", "Pearson.correlation.coefficient")]
matrix.plot.fre <- matrix.plot.fre %>%
  mutate(cell.type = factor(cell.type, levels = cell_order)) %>%
  arrange(cell.type, gene)

length(unique(matrix.plot.fre$cell.type))

matrix.plot.fre <- matrix.plot.fre %>%
  pivot_wider(names_from = cell.type,
              values_from = Pearson.correlation.coefficient) %>%
  as.data.frame()

rownames(matrix.plot.fre) <- matrix.plot.fre$gene
matrix.plot.fre <- matrix.plot.fre[ , -1]  


missing_cols <- setdiff(cell_order, colnames(matrix.plot.fre))
if(length(missing_cols) > 0){
  matrix.plot.fre[ , missing_cols] <- 0
}


matrix.plot.fre <- matrix.plot.fre[ , cell_order, drop = FALSE]
matrix.plot.fre[is.na(matrix.plot.fre)] <- 0

breaksList <- c(seq(0, 0.8, length.out=30)[-1],
                seq(0.8, 1, length.out=10)[-1])
my_palette <- colorRampPalette(c("white","#FCE8DD","#E9212C"))(length(breaksList)-1)
pheatmap(matrix.plot.fre,
         cluster_rows = F,
         cluster_cols = F,
         color = my_palette,
         border_color = NA,
         breaks = breaksList,
         main = paste0("6 cell types, ", gene.num," AIGs"), 
         filename = paste0("6.1.pheatmap.genes.ovelap.6 cell types, ", gene.num," AIGs.pdf"))


########################################################
######################down for plot#####################
########################################################
cor.matrix.down.data <- as.data.frame(table(cor.matrix.down$gene))
colnames(cor.matrix.down.data) <- c("Gene", "freq")

for (i in 1:5){
  gene.sel <- cor.matrix.down.data[cor.matrix.down.data$freq == i,]
  gene.num <- nrow(gene.sel)
  cell_order <- as.character(cor.plot.order$celltype)
  
  
  matrix.plot.fre <- cor.matrix.down[cor.matrix.down$gene %in% gene.sel$Gene,]
  matrix.plot.fre <- matrix.plot.fre[,c("gene",  "cell.type", "Pearson.correlation.coefficient")]
  matrix.plot.fre <- matrix.plot.fre %>%
    mutate(cell.type = factor(cell.type, levels = cell_order)) %>%
    arrange(cell.type, gene)
  
  length(unique(matrix.plot.fre$cell.type))
  
  matrix.plot.fre <- matrix.plot.fre %>%
    pivot_wider(names_from = cell.type,
                values_from = Pearson.correlation.coefficient) %>%
    as.data.frame()
  rownames(matrix.plot.fre) <- matrix.plot.fre$gene
  matrix.plot.fre <- matrix.plot.fre[ , -1]  # 去掉gene列
  
  missing_cols <- setdiff(cell_order, colnames(matrix.plot.fre))
  if(length(missing_cols) > 0){
    matrix.plot.fre[ , missing_cols] <- 0
  }
  

  matrix.plot.fre <- matrix.plot.fre[ , cell_order, drop = FALSE]
  matrix.plot.fre[is.na(matrix.plot.fre)] <- 0
  
  breaksList <- c(seq(-1, -0.8, length.out=10),
                  seq(-0.8, 0, length.out=30)[-1])
  my_palette <- colorRampPalette(c("#414986","#C2BDDB", "white" ))(length(breaksList)-1)
  
  pheatmap(matrix.plot.fre,
           cluster_rows = F,
           cluster_cols = F,
           color = my_palette,
           border_color = NA,
           breaks = breaksList,
           main = paste0(i," cell types, ", gene.num," AIGs"), 
           filename = paste0("6.1.pheatmap.genes.ovelap.",i," cell types, ", gene.num," ADGs.pdf"))
  
}

#################
nrow(cor.matrix.down.data[cor.matrix.down.data$freq > 5,] ) #284

###
gene.sel <- cor.matrix.down.data[cor.matrix.down.data$freq >5,]
write.table(gene.sel$Gene, file = "/share/home/yangjingyi/project/1.BD/up.load/1.normal/6.1.output/6.1.ADGs.6.cell.types.20250828.txt", quote = F,row.names = F,col.names = F)

gene.num <- nrow(gene.sel)

cell_order <- as.character(cor.plot.order$celltype)


matrix.plot.fre <- cor.matrix.down[cor.matrix.down$gene %in% gene.sel$Gene,]
matrix.plot.fre <- matrix.plot.fre[,c("gene",  "cell.type", "Pearson.correlation.coefficient")]
matrix.plot.fre <- matrix.plot.fre %>%
  mutate(cell.type = factor(cell.type, levels = cell_order)) %>%
  arrange(cell.type, gene)

length(unique(matrix.plot.fre$cell.type))

matrix.plot.fre <- matrix.plot.fre %>%
  pivot_wider(names_from = cell.type,
              values_from = Pearson.correlation.coefficient) %>%
  as.data.frame()
rownames(matrix.plot.fre) <- matrix.plot.fre$gene
matrix.plot.fre <- matrix.plot.fre[ , -1] 

missing_cols <- setdiff(cell_order, colnames(matrix.plot.fre))
if(length(missing_cols) > 0){
  matrix.plot.fre[ , missing_cols] <- 0
}

matrix.plot.fre <- matrix.plot.fre[ , cell_order, drop = FALSE]
matrix.plot.fre[is.na(matrix.plot.fre)] <- 0

breaksList <- c(seq(-1, -0.8, length.out=10),
                seq(-0.8, 0, length.out=30)[-1])
my_palette <- colorRampPalette(c("#414986","#C2BDDB", "white" ))(length(breaksList)-1)

pheatmap(matrix.plot.fre,
         cluster_rows = F,
         cluster_cols = F,
         color = my_palette,
         border_color = NA,
         breaks = breaksList,
         main = paste0("6 cell types, ", gene.num," ADGs"), 
         filename = paste0("6.1.pheatmap.genes.ovelap.6 cell types, ", gene.num," ADGs.pdf"))


#########################################################################################
####go.term##############################################################################
#########################################################################################
####
matrix <- read.table("/share/home/yangjingyi/project/1.BD/plot.final.v7/1.normal.cell.line.filt.organ/11.1.output/11.1.bp.AIGs.6.cell.types.20250828.txt", sep = "\t", header =T)
matrix <- matrix[matrix$PValue < 0.05,]
matrix <- matrix[order(matrix$PValue,decreasing = T),]
matrix$Term <- gsub("GO:[^~]*~","",matrix$Term)

matrix <- matrix[matrix$Term %in% c( "cell differentiation" ,                  
                                     "peptide cross-linking"  ,        
                                     "negative regulation of peptidase activity" ,
                                     "intermediate filament organization"  ,    
                                     "negative regulation of proteolysis" ),]
matrix$Term <- factor(matrix$Term, levels = matrix$Term)
p <-ggplot(matrix, aes(x = -log10(PValue), y = Term, fill = -log10(PValue))) +
  geom_col(width = 0.6) +
  scale_fill_gradient(low = "#FCE8DD", high = "#E9212C") +
  geom_text(
    aes(x = 0.05, label = Term),
    hjust = 0, size = 5
  ) +
  labs(x = expression(-log[10](p~value)), y = NULL,
       title = "Representative GO term enriched in AIGs") +
  theme(panel.background = element_rect(fill = 'white'), 
        axis.text.x = element_text(colour = "black",angle = 0, hjust = 1), 
        axis.text.y = element_blank(), 
        axis.line.x = element_line(colour = "black",  size=0.5, lineend = "butt"), 
        axis.line.y = element_line(colour = "black", size=0.5)) 
ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/plot.final.v7/1.normal.cell.line.filt.organ/11.1.output/11.1.bp.AIGs.6.cell.types.20250828.pdf",width = 8,height = 4)

####
matrix <- read.table("/share/home/yangjingyi/project/1.BD/plot.final.v7/1.normal.cell.line.filt.organ/11.1.output/11.1.bp.ADGs.6.cell.types.20250828.txt", sep = "\t", header =T)
matrix <- matrix[matrix$PValue < 0.05,]
matrix <- matrix[order(matrix$PValue,decreasing = T),]
matrix$Term <- gsub("GO:[^~]*~","",matrix$Term)

matrix <- matrix[matrix$Term %in% c( "telomere maintenance via telomerase" ,                  
                                     "hippo signaling"     ,        
                                     "cellular response to transforming growth factor beta stimulus"   ,
                                     "positive regulation of canonical Wnt signaling pathway",    
                                     "cell division" ),]
matrix$Term <- factor(matrix$Term, levels = matrix$Term)
p <-ggplot(matrix, aes(x = -log10(PValue), y = Term, fill = -log10(PValue))) +
  geom_col(width = 0.6) +
  scale_fill_gradient(low = "#C2BDDB", high = "#414986") +
  geom_text(
    aes(x = 0.05, label = Term),
    hjust = 0, size = 5
  ) +
  labs(x = expression(-log[10](p~value)), y = NULL,
       title = "Representative GO term enriched in ADGs") +
  theme(panel.background = element_rect(fill = 'white'), 
        axis.text.x = element_text(colour = "black",angle = 0, hjust = 1), 
        axis.text.y = element_blank(), 
        axis.line.x = element_line(colour = "black",  size=0.5, lineend = "butt"), 
        axis.line.y = element_line(colour = "black", size=0.5)) 
ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/plot.final.v7/1.normal.cell.line.filt.organ/11.1.output/11.1.bp.ADGs.6.cell.types.20250828.pdf",width = 8,height = 4)



########################################################################################
##########select genes##################################################################
########################################################################################
gene.sel <- cor.matrix.down.data[cor.matrix.down.data$freq == 10,]


gene.sel <- cor.matrix.up.data[cor.matrix.up.data$freq == 10,]

###
load("/share/home/yangjingyi/project/1.BD/up.load/1.normal/4.0.output/amniocytes.normal.Robj")
matrix <- as.matrix(GetAssayData(amniocytes.normal, slot = "data"))  
matrix <- as.data.frame(matrix[rownames(matrix) %in% c("SCP2", "PRSS27",  "TMPRSS2"), ])
matrix <- as.data.frame(t(matrix))
matrix$barcode <- rownames(matrix)

metadata <- amniocytes.normal@meta.data
metadata <- metadata[,c("barcode", "singleR_labels",  "orig.ident", "gestation.day")]

matrix <- merge(matrix,metadata, by="barcode")
df_plot <- matrix %>%
  group_by(gestation.day, singleR_labels) %>%
  summarise(mean_SCP2 = mean(SCP2, na.rm = TRUE),
            mean_PRSS27 = mean(PRSS27, na.rm = TRUE),
            mean_TMPRSS2 = mean(TMPRSS2, na.rm = TRUE),
            .groups = "drop") 
df_z <- df_plot %>%
  pivot_longer(cols = starts_with("mean_"),  
               names_to = "gene", 
               values_to = "expression") %>%
  group_by(singleR_labels, gene) %>%  
  mutate(zscore = scale(expression)[,1]) %>%  
  ungroup()

p <- ggplot(df_z[df_z$gene == "mean_SCP2",], aes(x = gestation.day, 
                                                 y = zscore, 
                                                 color = singleR_labels, 
                                                 group = singleR_labels)) +
  geom_vline(xintercept = seq(105, 217, by = 7),
             color = "grey85", linewidth = 0.4, linetype = "dashed",
             inherit.aes = FALSE) +
  geom_smooth(se = FALSE, method = "loess", size = 1) +
  labs(x = "Gestation day", y = "SCP2 expression", color = "Cell type")+
  scale_color_manual(values = c(
    "Intestine-Intestinal.epithelial.cells" ="#FC8DCA",
    "Kidney-Megakaryocytes" = "#8CCAFD",
    "Kidney-Ureteric.bud.cells" = "#009FF6",
    "Lung-Bronchiolar.and.alveolar.epithelial.cells" = "#FF662E",
    "Lung-Squamous.epithelial.cells" = "#FFC763",
    "Placenta-Megakaryocytes" = "#704D9E",
    "Stomach-Ciliated.epithelial.cells" = "#86CB66",
    "Stomach-Erythroblasts" = "#007A0B",
    "Stomach-Mesothelial.cells" = "#91AD5A",
    "Stomach-Myeloid.cells" = "#369891"
  )) +
  theme(axis.text.x = element_text(colour = "black", angle = 30, hjust = 1), 
        axis.text.y = element_text(colour = "black", hjust = 1),
        panel.background = element_rect(fill = 'white'), 
        axis.line.x = element_line(colour = "black",  size=0.5, lineend = "butt"), 
        axis.line.y = element_line(colour = "black", size=0.5)) 
ggsave(filename="6.1.ggplot.SCP2.expression.zscore.with.legend.20250529.pdf",p,width = 8, height = 8)


p <- ggplot(df_z[df_z$gene == "mean_PRSS27",], aes(x = gestation.day, 
                                                   y = zscore, 
                                                   color = singleR_labels, 
                                                   group = singleR_labels)) +
  geom_vline(xintercept = seq(105, 217, by = 7),
             color = "grey85", linewidth = 0.4, linetype = "dashed",
             inherit.aes = FALSE) +
  geom_smooth(se = FALSE, method = "loess", size = 1) +
  labs(x = "Gestation day", y = "PRSS27 expression", color = "Cell type")+
  scale_color_manual(values = c(
    "Intestine-Intestinal.epithelial.cells" ="#FC8DCA",
    "Kidney-Megakaryocytes" = "#8CCAFD",
    "Kidney-Ureteric.bud.cells" = "#009FF6",
    "Lung-Bronchiolar.and.alveolar.epithelial.cells" = "#FF662E",
    "Lung-Squamous.epithelial.cells" = "#FFC763",
    "Placenta-Megakaryocytes" = "#704D9E",
    "Stomach-Ciliated.epithelial.cells" = "#86CB66",
    "Stomach-Erythroblasts" = "#007A0B",
    "Stomach-Mesothelial.cells" = "#91AD5A",
    "Stomach-Myeloid.cells" = "#369891"
  )) +
  theme(axis.text.x = element_text(colour = "black", angle = 30, hjust = 1), 
        axis.text.y = element_text(colour = "black", hjust = 1),
        panel.background = element_rect(fill = 'white'), 
        axis.line.x = element_line(colour = "black",  size=0.5, lineend = "butt"), 
        axis.line.y = element_line(colour = "black", size=0.5)) 
ggsave(filename="6.1.ggplot.PRSS27.expression.zscore.with.legend.20250529.pdf",p,width = 8, height = 8)

p <- ggplot(df_z[df_z$gene == "mean_TMPRSS2",], aes(x = gestation.day, 
                                                    y = zscore, 
                                                    color = singleR_labels, 
                                                    group = singleR_labels)) +
  geom_vline(xintercept = seq(105, 217, by = 7),
             color = "grey85", linewidth = 0.4, linetype = "dashed",
             inherit.aes = FALSE) +
  geom_smooth(se = FALSE, method = "loess", size = 1) +
  labs(x = "Gestation day", y = "TMPRSS2 expression", color = "Cell type")+
  scale_color_manual(values = c(
    "Intestine-Intestinal.epithelial.cells" ="#FC8DCA",
    "Kidney-Megakaryocytes" = "#8CCAFD",
    "Kidney-Ureteric.bud.cells" = "#009FF6",
    "Lung-Bronchiolar.and.alveolar.epithelial.cells" = "#FF662E",
    "Lung-Squamous.epithelial.cells" = "#FFC763",
    "Placenta-Megakaryocytes" = "#704D9E",
    "Stomach-Ciliated.epithelial.cells" = "#86CB66",
    "Stomach-Erythroblasts" = "#007A0B",
    "Stomach-Mesothelial.cells" = "#91AD5A",
    "Stomach-Myeloid.cells" = "#369891"
  )) +
  theme(axis.text.x = element_text(colour = "black", angle = 30, hjust = 1), 
        axis.text.y = element_text(colour = "black", hjust = 1),
        panel.background = element_rect(fill = 'white'), 
        axis.line.x = element_line(colour = "black",  size=0.5, lineend = "butt"), 
        axis.line.y = element_line(colour = "black", size=0.5)) 
ggsave(filename="6.1.ggplot.TMPRSS2.expression.zscore.with.legend.20250529.pdf",p,width = 8, height = 8)
