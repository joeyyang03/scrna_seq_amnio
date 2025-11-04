rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/2.model/2.0.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib")

library(dplyr)
library(readr)
library(Seurat)
library(ggplot2)
library(pheatmap)
library(cowplot)
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/0.0.output/amnio.train.Robj")
folder_path <- "/share/home/yangjingyi/project/1.BD/up.load/2.model/1.0.output"
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

cor.matrix <-  data.frame(matrix(ncol = 7, nrow = 0)) 
colnames(cor.matrix) <- c("cell.type", "gene", "p.value", "Pearson.correlation.coefficient", "exp.max","exp.min", "exp.mean")
for (i in 1:length(csv_files)){
  res <- read.csv(csv_files[[i]],row.names = 1)
  cor.matrix <- rbind(cor.matrix,res)
}

save(cor.matrix, file = "/share/home/yangjingyi/project/1.BD/up.load/2.model/1.0.output/1.0.amniocytes.with.celltype.cor.exp.matrix.20251007.Robj")

dim(cor.matrix)

cor.matrix <- cor.matrix[cor.matrix$exp.mean > 0.1, ]
table(cor.matrix$cell.type)


##########################for plot
cor.plot.0.8 <- cor.matrix[abs(cor.matrix$Pearson.correlation.coefficient) > 0.8, ]
table(cor.plot.0.8$cell.type)

cor.plot.0.8$Correlation[cor.plot.0.8$Pearson.correlation.coefficient > 0.8 & cor.plot.0.8$p.value < 0.05] <- "Positive"
cor.plot.0.8$Correlation[cor.plot.0.8$Pearson.correlation.coefficient < -0.8 & cor.plot.0.8$p.value < 0.05] <- "Negative"
table(cor.plot.0.8$Correlation)

cor.plot <- as.data.frame(table(cor.plot.0.8$cell.type, cor.plot.0.8$Correlation))
colnames(cor.plot) <- c("celltype", "Correlation", "value")

cor.plot.order <- as.data.frame(table(cor.plot.0.8$cell.type))
colnames(cor.plot.order) <- c("celltype", "value")
cor.plot.order <- cor.plot.order[order(cor.plot.order$value, decreasing = T),]
cor.plot.order$celltype <- factor(cor.plot.order$celltype, levels = cor.plot.order$celltype)

cor.plot$celltype <- factor(cor.plot$celltype, levels = cor.plot.order$celltype)

p <- ggplot(cor.plot, aes(x = celltype, y = value, fill = Correlation)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Positive" = "brown", 
                               "Negative" = "steelblue"))+
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
    vjust = -0.3, hjust = 0.5)+ coord_flip() 
ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.0.output/2.0.cor.0.8.exp.mean.0.1.Numbers.of.genes.without.cut.20251007.pdf",width = 8,height = 8)

######################make df selection
cor.matrix.sel <- cor.matrix[cor.matrix$p.value < 0.05, ]
cor.matrix.sel <- cor.matrix.sel[abs(cor.matrix.sel$Pearson.correlation.coefficient) > 0.8, ]
table(cor.matrix.sel$cell.type)
cor.matrix.sel2 <- cor.matrix.sel %>%
  top_n(n = 30, 
        wt = abs(Pearson.correlation.coefficient))
table(cor.matrix.sel2$cell.type)

cor.gene <- unique(cor.matrix.sel2$gene)
cor.cell.type <- unique(cor.matrix.sel2$cell.type)


matrix <- as.data.frame(as.matrix(amniocytes.sample[["RNA"]]@data))
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

save(df, file = "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.0.output/2.0.amniocytes.with.celltype.cor.top30.exp.0.1.df.20251007.Robj")
#########plot
df.plot <- df

metadata <- amniocytes.sample@meta.data
metadata <- metadata[,c("gestation.day", "orig.ident")]
metadata <- metadata %>%
  group_by(orig.ident) %>%
  summarise(
    gestation.day = first(gestation.day)
  )
metadata <- as.data.frame(metadata)
df.plot <- merge(df.plot,metadata,by = "orig.ident")
df.plot[is.na(df.plot)] <- 0

df.plot <- df.plot[order(df.plot$gestation.day,decreasing = F),]
df.plot <- df.plot %>%dplyr::select(-orig.ident, -gestation.day)

df.plot <- as.data.frame(t(df.plot))

my_colors <- colorRampPalette(c("steelblue", "white", "brown"))(100)


my_breaks <- seq(-2, 2, length.out = 101)  
pheatmap(df.plot,
         scale = "row",
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         border_color = NA,
         show_colnames = FALSE,
         color = my_colors,
         breaks = my_breaks,
         filename = "2.0.pheatmap.df.train.20251007.pdf")
####plot
cor.celltype <- unique(cor.matrix$cell.type)
range(-log10(cor.matrix$p.value))


global_p_range <- range(-log10(cor.matrix$p.value))
global_size_breaks <- pretty(global_p_range, n = 5) 
p <- list()
for (i in 1:length(cor.celltype)){
  celltype <- cor.celltype[[i]]
  res <- cor.matrix[cor.matrix$cell.type == celltype,]
  res <- res[order(res$Pearson.correlation.coefficient, decreasing = F),]
  res$Correlation[res$Pearson.correlation.coefficient > 0.8 & res$p.value < 0.05] <- "Positive"
  res$Correlation[res$Pearson.correlation.coefficient < -0.8 & res$p.value < 0.05] <- "Negative"
  res$Correlation[abs(res$Pearson.correlation.coefficient) < 0.8] <- "Uncorrelated"
  res$gene <- factor(res$gene,levels =res$gene)
  
  
  p[[i]] <- ggplot(res, 
                   aes(x = gene, y = Pearson.correlation.coefficient,
                       size = -log10(p.value), color = Correlation)) +
    geom_point(alpha = 1) +
    scale_color_manual(values = c("Positive" = "brown", 
                                  "Uncorrelated" = "grey",
                                  "Negative" = "steelblue")) +
    scale_size_continuous(
      range = c(0.02, 4),   
      limits = c(0, 16),   
      breaks = c(1,2,10,15),
      guide = guide_legend(
        override.aes = list(shape = 16, alpha = 1),
        nrow = 1
      )
    ) +
    coord_cartesian(clip = "off") +
    scale_x_discrete(expand = expansion(mult = 0.03)) +
    ylim(-1, 1) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.background = element_blank(), 
      axis.line.x = element_line(colour = "black", size = 0.5), 
      axis.line.y = element_line(colour = "black", size = 0.5)
    ) +
    labs(title = celltype, x = "gene", y = "Correlation.coefficient")
}
plot <- plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],
                  p[[7]],p[[8]],p[[9]],p[[10]],ncol = 3)
ggsave(filename="2.1.celltype.cor.0.8.exp.mean.0.1.with.legend.20250711.pdf",plot,width = 10, height = 9)

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
save(df.test, file = "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.0.output/2.0.amniocytes.test.with.celltype.cor.top30.exp.0.1.df.20251007.Robj")
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
save(df.FGR.cohort1, file = "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.0.output/2.0.amniocytes.FGR.cohort1.with.celltype.cor.top30.exp.0.1.df.20251007.Robj")
