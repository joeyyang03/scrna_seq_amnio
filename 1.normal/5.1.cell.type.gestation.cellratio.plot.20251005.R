rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/1.normal/5.1.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 

library(Seurat)
library(ggplot2)
library(tidyverse)
library(rhdf5)
library("edgeR")
library(colorspace)
library(RColorBrewer)
library(data.table)
library(dplyr)
library(cowplot)
library(ggpubr)
library(patchwork)
library(pheatmap)
load("/share/home/yangjingyi/project/1.BD/up.load/1.normal/4.0.output/amniocytes.normal.Robj")
###
metadata <- amniocytes.normal@meta.data
metadata$gestation.week <- NA
metadata$gestation.week[which(str_detect(metadata$orig.ident, "^GW15"))] <- "GW15"
metadata$gestation.week[which(str_detect(metadata$orig.ident, "^GW17"))] <- "GW17"
metadata$gestation.week[which(str_detect(metadata$orig.ident, "^GW18"))] <- "GW18"
metadata$gestation.week[which(str_detect(metadata$orig.ident, "^GW19"))] <- "GW19"
metadata$gestation.week[which(str_detect(metadata$orig.ident, "^GW20"))] <- "GW20"
metadata$gestation.week[which(str_detect(metadata$orig.ident, "^GW21"))] <- "GW21"
metadata$gestation.week[which(str_detect(metadata$orig.ident, "^GW22"))] <- "GW22"
metadata$gestation.week[which(str_detect(metadata$orig.ident, "^GW23"))] <- "GW23"
metadata$gestation.week[which(str_detect(metadata$orig.ident, "^GW24"))] <- "GW24"
metadata$gestation.week[which(str_detect(metadata$orig.ident, "^GW25"))] <- "GW25"
metadata$gestation.week[which(str_detect(metadata$orig.ident, "^GW26"))] <- "GW26"
metadata$gestation.week[which(str_detect(metadata$orig.ident, "^GW27"))] <- "GW27"
metadata$gestation.week[which(str_detect(metadata$orig.ident, "^GW28"))] <- "GW28"
metadata$gestation.week[which(str_detect(metadata$orig.ident, "^GW29"))] <- "GW29"
metadata$gestation.week[which(str_detect(metadata$orig.ident, "^GW30"))] <- "GW30"
table(metadata$gestation.week,metadata$orig.ident)
amniocytes.normal@meta.data <- metadata 

############
###
metadata$organ <- NA
metadata$organ[which(str_detect(metadata$singleR_labels, "^Stomach"))] <- "Stomach"
metadata$organ[which(str_detect(metadata$singleR_labels, "^Kidney"))] <- "Kidney"
metadata$organ[which(str_detect(metadata$singleR_labels, "^Intestine"))] <- "Intestine"
metadata$organ[which(str_detect(metadata$singleR_labels, "^Lung"))] <- "Lung"
metadata$organ[which(str_detect(metadata$singleR_labels, "^Placenta"))] <- "Placenta"
amniocytes.normal@meta.data <- metadata

###
amniocytes.gestation.cell.type <- table(metadata$organ,metadata$gestation.week)
amniocytes.gestation.cell.type <- matrix(amniocytes.gestation.cell.type, 
                                         ncol=ncol(amniocytes.gestation.cell.type), 
                                         dimnames=dimnames(amniocytes.gestation.cell.type))

amniocytes.gestation.cell.type <- as.data.frame(amniocytes.gestation.cell.type)
amniocytes.gestation.cell.type <- sweep(amniocytes.gestation.cell.type, 2, colSums(amniocytes.gestation.cell.type), `/`)
apply(amniocytes.gestation.cell.type, 2, sum)


pheatmap(amniocytes.gestation.cell.type,
         scale = "row",color = colorRampPalette(c("#a4cde1", "white", "#af93c4" ))(100),
         cluster_rows = T,cluster_cols = F,
         border_color = F, filename = "5.1.pheatmap.gestation.organ.ratio.20251005.pdf")

#############################
amniocytes.gestation.cell.type <- table(metadata$singleR_labels, metadata$gestation.week)
amniocytes.gestation.cell.type <- matrix(amniocytes.gestation.cell.type, 
                                         ncol=ncol(amniocytes.gestation.cell.type), 
                                         dimnames=dimnames(amniocytes.gestation.cell.type))

amniocytes.gestation.cell.type <- as.data.frame(amniocytes.gestation.cell.type)
amniocytes.gestation.cell.type <- sweep(amniocytes.gestation.cell.type, 2, colSums(amniocytes.gestation.cell.type), `/`)
amniocytes.gestation.cell.type$cell.type <- rownames(amniocytes.gestation.cell.type)
amniocytes.gestation.cell.type <- reshape2::melt(amniocytes.gestation.cell.type, id="cell.type")

#####
df <- amniocytes.gestation.cell.type

df$Gestation.week <- as.numeric(gsub("GW", "", df$variable))
res_list <- lapply(split(df, df$cell.type), function(d) {
  fit <- lm(value ~ Gestation.week, data = d)
  pval <- summary(fit)$coefficients["Gestation.week", "Pr(>|t|)"]
  data.frame(cell.type = unique(d$cell.type), pval = pval)
})

res <- do.call(rbind, res_list)

# FDR 
res$FDR <- p.adjust(res$pval, method = "BH")
res$logFDR <- -log10(res$FDR)
res$Sig <- res$FDR < 0.05
head(res)
res <- res[order(res$FDR, decreasing = F),]
res$cell.type <- factor(res$cell.type, levels = res$cell.type)

p1 <- ggplot(res, aes(x = cell.type, y = logFDR, fill = Sig)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE"= "#af93c4", "FALSE"= "#a4cde1")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


amniocytes.gestation.cell.type$Sig <- NA
amniocytes.gestation.cell.type$Sig[amniocytes.gestation.cell.type$cell.type %in% c("Intestine-Intestinal.epithelial.cells",
                                                                                   "Kidney-Ureteric.bud.cells",
                                                                                   "Placenta-Megakaryocytes",
                                                                                   "Stomach-Mesothelial.cells"  ,
                                                                                   "Stomach-Myeloid.cells")] <- "TRUE"
amniocytes.gestation.cell.type$Sig[amniocytes.gestation.cell.type$cell.type %in% c("Kidney-Megakaryocytes",
                                                                                   "Lung-Bronchiolar.and.alveolar.epithelial.cells",
                                                                                   "Lung-Squamous.epithelial.cells",
                                                                                   "Stomach-Ciliated.epithelial.cells"  ,
                                                                                   "Stomach-Erythroblasts")] <- "FALSE"
amniocytes.gestation.cell.type$cell.type <- factor(amniocytes.gestation.cell.type$cell.type, levels = res$cell.type)

p2 <- ggplot(amniocytes.gestation.cell.type, 
             aes(x = cell.type, y = variable, 
                 size = value*100, color = Sig)) +  
  geom_point() +
  scale_size_continuous(range = c(0.5, 10), breaks = c(10, 20, 40)) +
  scale_color_manual(values = c("TRUE"= "#af93c4", "FALSE"= "#a4cde1")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

final_plot <- p1 / p2 + plot_layout(heights = c(1, 3))
final_plot
ggsave(filename="5.1.fdr.celltype.ratio.20251005.pdf",final_plot,width = 6, height = 10)


###############
celltype.sel <- res[res$FDR < 0.05,]$cell.type

amniocytes.celltype.ratio <- table(metadata$singleR_labels, metadata$orig.ident)
amniocytes.celltype.ratio <- matrix(amniocytes.celltype.ratio, 
                                    ncol=ncol(amniocytes.celltype.ratio), 
                                    dimnames=dimnames(amniocytes.celltype.ratio))

amniocytes.celltype.ratio <- as.data.frame(amniocytes.celltype.ratio)
amniocytes.celltype.ratio <- sweep(amniocytes.celltype.ratio, 2, colSums(amniocytes.celltype.ratio), `/`)
amniocytes.celltype.ratio$cell.type <- rownames(amniocytes.celltype.ratio)
amniocytes.celltype.ratio <- reshape2::melt(amniocytes.celltype.ratio, id="cell.type")

celltype.sel <- unique(res$cell.type)
p <- list()
for (i in 1:length(celltype.sel)){
  a <- celltype.sel[i]
  data <- amniocytes.celltype.ratio[amniocytes.celltype.ratio$cell.type == a,]
  data$Gestation.week <- as.numeric(gsub("GW(\\d+)_.*", "\\1", data$variable))	
  data$Gestation.week.factor <- factor(data$Gestation.week,
                                       levels = sort(unique(data$Gestation.week)))
  data$GestationNum <- data$Gestation.week
  

  p[[i]] <- ggplot(data, aes(x = Gestation.week.factor, y = value)) +
    geom_boxplot(aes(fill = Gestation.week.factor), fill = "#704D9E", width = 0.5, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.8) +
    geom_smooth(data = data, aes(x = as.numeric(Gestation.week.factor), y = value),
                method = "lm", color = "black", se = TRUE, inherit.aes = FALSE) +
    stat_cor(data = data, aes(x = GestationNum, y = value),
             method = "pearson", label.x = 5, label.y = max(data$value, na.rm = TRUE)) +
    theme(
      panel.background = element_rect(fill = "white", color = "black", linewidth = 0.6), 
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),   
      axis.line = element_blank(),    
      axis.text.x = element_text(color = "black", angle = 30, hjust = 1, size = 10),
      axis.text.y = element_text(color = "black", size = 10),
      axis.title = element_text(face = "bold", size = 11),
      plot.title = element_text(hjust = 0, face = "bold", size = 12),
      legend.position = "none",
      plot.margin = margin(6, 6, 6, 6)
    )+ggtitle(a)
}



plot <- plot_grid(p[[1]],
                  p[[2]],
                  p[[4]],
                  p[[3]],
                  p[[5]],
                  ncol=3)
plot
ggsave(filename="5.1.celltype.ratio.lm.with.pvalue.significant.20251005.pdf",plot,width = 10, height = 5)

plot <- plot_grid(p[[6]],
                  p[[7]],
                  p[[9]],
                  p[[8]],
                  p[[10]],
                  ncol=3)
plot
ggsave(filename="5.1.celltype.ratio.lm.without.pvalue.significant.20251005.pdf",plot,width = 10, height = 5)
