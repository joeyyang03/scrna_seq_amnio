rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/1.normal/7.2.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 

library(dplyr)
library("Mfuzz")
library(Seurat)
library(tibble)
library(Cairo)
library(ggplot2)
library(pheatmap)
library(ggridges)
folder_path <- "/share/home/yangjingyi/project/1.BD/up.load/1.normal/7.1.output"
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

wave.matrix <-  data.frame(matrix(ncol = 5, nrow = 0)) 
colnames(wave.matrix) <- c("Cluster",  "Gene", "Max_Expression", "Peak_Time", "cell.type" )
for (i in 1:length(csv_files)){
  res <- read.csv(csv_files[[i]],row.names = 1)
  wave.matrix <- rbind(wave.matrix,res)
  
}
######################################wave.summary
wave.summary.plot <- as.data.frame(table(wave.matrix$Peak_Time))
colnames(wave.summary.plot) <- c("Peak_Time", "Count")

new_row <- data.frame(Peak_Time = "GW23", Count = 0)

wave.summary.plot <- rbind(wave.summary.plot, new_row)
wave.summary.plot$Week <- as.numeric(gsub("GW", "", wave.summary.plot$Peak_Time))

wave.summary.plot$Peak_Time <- factor(wave.summary.plot$Peak_Time, levels = paste0("GW", sort(unique(wave.summary.plot$Week))))
wave.summary.plot <- wave.matrix[,c("Peak_Time"),drop=F]
wave.summary.plot$Week <- as.numeric(gsub("GW", "", wave.summary.plot$Peak_Time))
p <- ggplot(wave.summary.plot, aes(x = Week)) +
  geom_density(fill = "#AF93C4", alpha = 0.6, color = NA) +
  scale_x_continuous(breaks = c(15, 17, 19, 21, 23, 25, 27, 29)) +
  labs(
    title = "Peak Time Density",
    x = "Gestation Week",
    y = "Density"
  ) +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),    
    axis.line = element_blank(),     
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(2, "mm"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 11),
    plot.title = element_text(
      hjust = 0, face = "bold", size = 13, margin = margin(b = 8)
    ),
    plot.margin = margin(8, 8, 8, 8)
  )
ggsave(p,file= "/share/home/yangjingyi/project/1.BD/up.load/1.normal/7.2.output/7.2.densityplot.total.gene.num.gestation.20251007.pdf",width = 8,height = 8)

#########################Ridge plot
table(wave.matrix$cell.type,wave.matrix$Peak_Time)

####
gestation.cell.type <- table(wave.matrix$cell.type,wave.matrix$Peak_Time)
gestation.cell.type <- matrix(gestation.cell.type, 
                              ncol=ncol(gestation.cell.type), 
                              dimnames=dimnames(gestation.cell.type))
gestation.cell.type <- as.data.frame(gestation.cell.type)
gestation.cell.type$GW23 <- 0
gestation.cell.type$celltype <- rownames(gestation.cell.type)

gestation.cell.type <- reshape2::melt(gestation.cell.type, id="celltype")

cell.type <- unique(gestation.cell.type$celltype)
max(gestation.cell.type$value)


gestation.cell.type$variable <- factor(
  gestation.cell.type$variable,
  levels = c("GW15","GW17","GW18","GW19","GW20","GW21",
             "GW22","GW23","GW24","GW25","GW26","GW27",
             "GW28","GW29","GW30")
)
#########################################################################################################################
wave.plot <- wave.matrix[,c("cell.type", "Peak_Time")]
wave.plot$Peak_Time <-  as.numeric(gsub("GW", "", as.character(wave.plot$Peak_Time)))
wave.plot$cell.type <- factor(wave.plot$cell.type,
                levels=c("Placenta-Megakaryocytes","Intestine-Intestinal.epithelial.cells" ,  "Lung-Bronchiolar.and.alveolar.epithelial.cells",
               "Lung-Squamous.epithelial.cells" ,  "Kidney-Ureteric.bud.cells", "Kidney-Megakaryocytes",                    
               "Stomach-Ciliated.epithelial.cells",   "Stomach-Mesothelial.cells" ,   "Stomach-Erythroblasts",  
               "Stomach-Myeloid.cells"   ))

p <- ggplot(wave.plot, 
            aes(x = Peak_Time, 
                y = cell.type, 
                fill = ..x..)) +
  geom_density_ridges_gradient(
    alpha = 0.7,
    scale = 0.9,
    rel_min_height = 0.01,
    bandwidth = 1,
    color = "black",      
    baseline_col = "black",
    size = 0.3
  ) +
  scale_x_continuous(breaks = seq(15, 30, by = 2)) +
  scale_fill_gradientn(
    colours = c("#F5A551", "#A277AF"),
    name = "Gest Week",
    guide = guide_colorbar(reverse = TRUE)
  ) +
  labs(
    x = "Gestation Week",
    y = NULL,
    title = "Cell Type Temporal Distribution"
  ) +
  theme_ridges(font_size = 11) +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 5)),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0, face = "bold", size = 13, margin = margin(b = 8)),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(8, 8, 8, 8)
  )
ggsave(p,file= "/share/home/yangjingyi/project/1.BD/up.load/1.normal/7.2.output/7.2.ridgeplot.celltype.gene.num.gestation.20251007.pdf",width = 6,height = 8)



early.mid <- c("GW15", "GW17")
mid <- c("GW18", "GW19", "GW20", "GW21", "GW22", "GW23", "GW24", "GW25", "GW26", "GW27")
late <- c("GW28", "GW29", "GW30")

celltype <- unique(wave.matrix$cell.type)
for (i in 1:length(celltype)){
  a <- celltype[[i]]
  res <- wave.matrix[wave.matrix$cell.type == a,]
  res1 <- res[res$Peak_Time %in% early.mid,]
  gene <- res1$Gene
  gene <- gsub("Mean_","",gene)
  gene <- unique(gene)
  write.table(gene, file = paste0("/share/home/yangjingyi/project/1.BD/up.load/1.normal/7.2.output/7.2.",a,".early.mid.20251007.txt"), quote = F,row.names = F,col.names = F)
  
  res2 <- res[res$Peak_Time %in% mid,]
  gene <- res2$Gene
  gene <- unique(gsub("Mean_","",gene))
  write.table(gene, file = paste0("/share/home/yangjingyi/project/1.BD/up.load/1.normal/7.2.output/7.2.",a,".mid.20251007.txt"), quote = F,row.names = F,col.names = F)
  
  res3 <- res[res$Peak_Time %in% late,]
  gene <- res3$Gene
  gene <- unique(gsub("Mean_","",gene))
  write.table(gene, file = paste0("/share/home/yangjingyi/project/1.BD/up.load/1.normal/7.2.output/7.2.",a,".late.20251007.txt"), quote = F,row.names = F,col.names = F)
  
}

