rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/5.1.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 

library(Seurat)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(reshape2)
library(ggpubr)
"%!in%" <- function(x,y)!('%in%'(x,y))
gestation.day <- read.csv("/share/home/yangjingyi/project/1.BD/data/rawdata.20251001/df.summary.FGR.20251008.csv")
colnames(gestation.day) <- c("Maternal.age",   "Gestation.day",  "indicator",      "Frozen",         "vali.group",     "gestation.week", "orig.ident" )

gestation.day.cohort1 <- gestation.day[gestation.day$vali.group == "cohort1",]
class(gestation.day.cohort1$Gestation.day)
gestation.day.cohort1 <- gestation.day.cohort1[order(gestation.day.cohort1$Gestation.day, decreasing = T),]
gestation.day.cohort1$orig.ident <- factor(gestation.day.cohort1$orig.ident, levels = gestation.day.cohort1$orig.ident)
p <- ggplot(gestation.day.cohort1, aes(x = orig.ident, y = Gestation.day)) +
  geom_hline(yintercept = seq(154, 232, by = 7),
             color = "grey85", linewidth = 0.4, linetype = "dashed",
             inherit.aes = FALSE) +
  geom_segment(aes(xend = orig.ident, y = min(Gestation.day) -16, yend = Gestation.day), color = "#4F4DA5") +
  geom_point(size = 1, shape = 19, stroke = 1.5, color = "#4F4DA5") +
  labs(title = "Gestation day (cohort1)",
       x = "sample",
       y = "Gestation.day") +
  coord_flip()+
  scale_y_continuous(expand = c(0,0), limits=c(140,240))+
  theme(
    panel.background = element_rect(fill = 'white'),
    axis.text.y = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black"),
    axis.line.x = element_line(colour = "black", size = 0.5),
    axis.line.y = element_line(colour = "black", size = 0.5)
  )
ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/up.load/3.FGR/5.1.output/5.1.lollipop.gestation.day.cohort1.20251007.pdf",width = 8,height = 8)

###
gestation.day.cohort2 <- gestation.day[gestation.day$vali.group == "cohort2",]
gestation.day.cohort2 <- gestation.day.cohort2[order(gestation.day.cohort2$Gestation.day, decreasing = T),]
gestation.day.cohort2$orig.ident <- factor(gestation.day.cohort2$orig.ident, levels = gestation.day.cohort2$orig.ident)
p <- ggplot(gestation.day.cohort2, aes(x = orig.ident, y = Gestation.day)) +
  geom_hline(yintercept = seq(126, 225, by = 7),
             color = "grey85", linewidth = 0.4, linetype = "dashed",
             inherit.aes = FALSE) +
  geom_segment(aes(xend = orig.ident, y = min(Gestation.day) -7, yend = Gestation.day), color = "#B366C7") +
  geom_point(size = 1, shape = 19, stroke = 1.5, color = "#B366C7") +
  labs(title = "Gestation day (cohort2)",
       x = "sample",
       y = "Gestation.day") +
  coord_flip()+
  scale_y_continuous(expand = c(0,0), limits=c(120,230))+
  theme(
    panel.background = element_rect(fill = 'white'),
    axis.text.y = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black"),
    axis.line.x = element_line(colour = "black", size = 0.5),
    axis.line.y = element_line(colour = "black", size = 0.5)
  )
ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/up.load/3.FGR/5.1.output/5.1.lollipop.gestation.day.cohort2.20251007.pdf",width = 8,height = 8)

load("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/5.0.output/amniocytes.cohort1.FGR.Robj")
amniocytes.FGR.cohort1@meta.data$singleR_labels <- factor(amniocytes.FGR.cohort1@meta.data$singleR_labels, levels = c("Stomach-Myeloid.cells",                         
                                                                                                                      "Stomach-Erythroblasts",                         
                                                                                                                      "Stomach-Mesothelial.cells",                     
                                                                                                                      "Stomach-Ciliated.epithelial.cells", 
                                                                                                                      "Kidney-Megakaryocytes", 
                                                                                                                      "Kidney-Ureteric.bud.cells", 
                                                                                                                      "Lung-Squamous.epithelial.cells",                
                                                                                                                      "Lung-Bronchiolar.and.alveolar.epithelial.cells",
                                                                                                                      "Placenta-Megakaryocytes",                       
                                                                                                                      "Intestine-Intestinal.epithelial.cells"))

color.cell.type <- c("#369891", "#86CB66","#91AD5A", "#007A0B", "#8CCAFD", "#009FF6", "#FF662E", "#FFC763", "#704D9E", "#FC8DCA")

p <- DimPlot(object = amniocytes.FGR.cohort1, reduction = "umap", group.by = "singleR_labels", label = F, repel = TRUE, cols = color.cell.type)
ggsave (p,file= "/share/home/yangjingyi/project/1.BD/up.load/3.FGR/5.1.output/5.1.dimplot.amniocytes.cohort1.FGR.celltype.with.legend.20251007.pdf",width = 8,height = 8)

#####
load("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/5.0.output/amniocytes.cohort2.FGR.Robj")
amniocytes.FGR.cohort2@meta.data$singleR_labels <- factor(amniocytes.FGR.cohort2@meta.data$singleR_labels, levels = c("Stomach-Myeloid.cells",                         
                                                                                                                      "Stomach-Erythroblasts",                         
                                                                                                                      "Stomach-Mesothelial.cells",                     
                                                                                                                      "Stomach-Ciliated.epithelial.cells", 
                                                                                                                      "Kidney-Megakaryocytes", 
                                                                                                                      "Kidney-Ureteric.bud.cells", 
                                                                                                                      "Lung-Squamous.epithelial.cells",                
                                                                                                                      "Lung-Bronchiolar.and.alveolar.epithelial.cells",
                                                                                                                      "Placenta-Megakaryocytes",                       
                                                                                                                      "Intestine-Intestinal.epithelial.cells"))

color.cell.type <- c("#369891", "#86CB66","#91AD5A", "#007A0B", "#8CCAFD", "#009FF6", "#FF662E", "#FFC763", "#704D9E", "#FC8DCA")

p <- DimPlot(object = amniocytes.FGR.cohort2, reduction = "umap", group.by = "singleR_labels", label = F, repel = TRUE, cols = color.cell.type)
ggsave (p,file= "/share/home/yangjingyi/project/1.BD/up.load/3.FGR/5.1.output/5.1.dimplot.amniocytes.cohort2.FGR.celltype.with.legend.20251007.pdf",width = 8,height = 8)
