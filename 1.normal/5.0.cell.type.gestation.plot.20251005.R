rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/1.normal/5.0.output/")
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
load("/share/home/yangjingyi/project/1.BD/up.load/1.normal/4.0.output/amniocytes.normal.Robj")

gestation.day <- read.csv("/share/home/yangjingyi/project/1.BD/data/rawdata.20251001/df.summary.normal.20251001.csv")
gestation.day <- gestation.day[gestation.day$Group == "Normal",]
class(gestation.day$Gestation.day)
colnames(gestation.day) <- c("Gestation.day"  ,           "orig.ident"   ,              "Group"  )
gestation.day <- gestation.day[order(gestation.day$Gestation.day, decreasing = T),]
gestation.day$orig.ident <- factor(gestation.day$orig.ident, levels = gestation.day$orig.ident)
max(gestation.day$Gestation.day)
min(gestation.day$Gestation.day)
p <- ggplot(gestation.day, aes(x = orig.ident, y = Gestation.day)) +
  geom_hline(yintercept = seq(105, 217, by = 7),
             color = "grey85", linewidth = 0.4, linetype = "dashed",
             inherit.aes = FALSE) +
  geom_segment(aes(xend = orig.ident, y = min(Gestation.day)-7, yend = Gestation.day), color = "#F5A551") + 
  geom_point(size = 1, shape = 19, stroke = 1.5, color = "#F5A551") +
  labs(title = "Gestation day (Normal.group)",
       x = "sample",
       y = "Gestation.day") +
  coord_flip()+
  scale_y_continuous(expand = c(0,0), limits=c(100,230))+
  theme(
    panel.background = element_rect(fill = 'white'),
    axis.text.y = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black"),
    axis.line.x = element_line(colour = "black", size = 0.5),
    axis.line.y = element_line(colour = "black", size = 0.5)
  )
ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/up.load/1.normal/5.0.output/5.0.lollipop.gestation.day.normal.group.20251005.pdf",width = 8,height = 8)


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

#################################################
###Gestation.age################################
###############################################
table(metadata$gestation.week)
amniocytes.sample <- table(metadata$gestation.week)
amniocytes.sample <- as.data.frame(amniocytes.sample)
colnames(amniocytes.sample) <- c("gestation.age", "Freq")

p <- ggplot(data=amniocytes.sample)+
  geom_bar(stat = 'identity', mapping=aes(x = gestation.age, y = Freq), fill = "#F5A551")+
  scale_y_continuous(expand = c(0,0), limits=c(0,58000))+
  theme(axis.text.x = element_text(colour = "black", angle = 30, hjust = 1), 
        axis.text.y = element_text(colour = "black", hjust = 1),
        legend.position = "none", panel.background = element_rect(fill = 'white'), 
        axis.line.x = element_line(colour = "black",  size=0.5, lineend = "butt"), 
        axis.line.y = element_line(colour = "black", size=0.5)) 
ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/up.load/1.normal/5.0.output/5.0.ggplot.gestation.week.normal.group.20251005.pdf",width = 8,height = 8)

###################################
######Cell.type#######################
####################################
unique(amniocytes.normal@meta.data$singleR_labels)
 
#########cell.line.ratio
metadata$group <- NA
metadata$group[which(metadata$singleR_labels %in% c("Stomach-Mesothelial.cells",                     
                                                    "Stomach-Ciliated.epithelial.cells", 
                                                    "Kidney-Ureteric.bud.cells", 
                                                    "Lung-Squamous.epithelial.cells",                
                                                    "Lung-Bronchiolar.and.alveolar.epithelial.cells",
                                                    "Intestine-Intestinal.epithelial.cells"))] <- "Epithelial_cells"
metadata$group[which(metadata$singleR_labels %in% c("Stomach-Myeloid.cells",                         
                                                    "Kidney-Megakaryocytes", 
                                                    "Placenta-Megakaryocytes",                       
                                                    "Stomach-Erythroblasts"))] <- "Hematopoietic_cells"

table(metadata$group, metadata$gestation.week)

amniocytes.gestation.cell.line <- table(metadata$group, metadata$gestation.week)
amniocytes.gestation.cell.line <- matrix(amniocytes.gestation.cell.line, 
                                         ncol=ncol(amniocytes.gestation.cell.line), 
                                         dimnames=dimnames(amniocytes.gestation.cell.line))
amniocytes.gestation.cell.line <- as.data.frame(amniocytes.gestation.cell.line)
amniocytes.gestation.cell.line <- sweep(amniocytes.gestation.cell.line, 2, colSums(amniocytes.gestation.cell.line), `/`)
amniocytes.gestation.cell.line$cell.line <- rownames(amniocytes.gestation.cell.line)

amniocytes.gestation.cell.line <- reshape2::melt(amniocytes.gestation.cell.line, id="cell.line")
amniocytes.gestation.cell.line$variable <- factor(amniocytes.gestation.cell.line$variable, levels = c("GW30", 
                                                                                                      "GW29",
                                                                                                      "GW28",
                                                                                                      "GW27", 
                                                                                                      "GW26", 
                                                                                                      "GW25", 
                                                                                                      "GW24", 
                                                                                                      "GW23", 
                                                                                                      "GW22", 
                                                                                                      "GW21", 
                                                                                                      "GW20", 
                                                                                                      "GW19",
                                                                                                      "GW18", 
                                                                                                      "GW17",
                                                                                                      "GW15"))

p <- ggplot(data=amniocytes.gestation.cell.line, mapping=aes(y=variable, x=value, fill=cell.line))+
  geom_bar(stat = "identity", position = "fill")+scale_fill_manual(values = c("#b2db87","#a48cbe"))+ 
  scale_x_continuous(expand = c(0,0))+
  theme(panel.background = element_rect(fill = 'white'), 
        axis.text.x = element_text(colour = "black",angle = 0, hjust = 1), 
        axis.text.y = element_text(colour = "black",angle = 0, hjust = 1), 
        axis.line.x = element_line(colour = "black",  size=0.5, lineend = "butt"), 
        axis.line.y = element_line(colour = "black", size=0.5)) 
ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/up.load/1.normal/5.0.output/5.0.amniocytes.cellline.gestation.ratio.20251005.pdf",width = 8,height = 8)
####

amniocytes.gestation.cell.line <- table(metadata$group, metadata$orig.ident)
amniocytes.gestation.cell.line <- matrix(amniocytes.gestation.cell.line, 
                                         ncol=ncol(amniocytes.gestation.cell.line), 
                                         dimnames=dimnames(amniocytes.gestation.cell.line))
amniocytes.gestation.cell.line <- as.data.frame(amniocytes.gestation.cell.line)
amniocytes.gestation.cell.line <- sweep(amniocytes.gestation.cell.line, 2, colSums(amniocytes.gestation.cell.line), `/`)
amniocytes.gestation.cell.line$cell.line <- rownames(amniocytes.gestation.cell.line)

amniocytes.gestation.cell.line <- reshape2::melt(amniocytes.gestation.cell.line, id="cell.line")
amniocytes.gestation.cell.line$stage <- NA
amniocytes.gestation.cell.line$stage[which(str_detect(amniocytes.gestation.cell.line$variable, "^GW15"))] <- "Stage I (GW15-GW27)"
amniocytes.gestation.cell.line$stage[which(str_detect(amniocytes.gestation.cell.line$variable, "^GW17"))] <- "Stage I (GW15-GW27)"
amniocytes.gestation.cell.line$stage[which(str_detect(amniocytes.gestation.cell.line$variable, "^GW18"))] <- "Stage I (GW15-GW27)"
amniocytes.gestation.cell.line$stage[which(str_detect(amniocytes.gestation.cell.line$variable, "^GW19"))] <- "Stage I (GW15-GW27)"
amniocytes.gestation.cell.line$stage[which(str_detect(amniocytes.gestation.cell.line$variable, "^GW20"))] <- "Stage I (GW15-GW27)"
amniocytes.gestation.cell.line$stage[which(str_detect(amniocytes.gestation.cell.line$variable, "^GW21"))] <- "Stage I (GW15-GW27)"
amniocytes.gestation.cell.line$stage[which(str_detect(amniocytes.gestation.cell.line$variable, "^GW22"))] <- "Stage I (GW15-GW27)"
amniocytes.gestation.cell.line$stage[which(str_detect(amniocytes.gestation.cell.line$variable, "^GW23"))] <- "Stage I (GW15-GW27)"
amniocytes.gestation.cell.line$stage[which(str_detect(amniocytes.gestation.cell.line$variable, "^GW24"))] <- "Stage I (GW15-GW27)"
amniocytes.gestation.cell.line$stage[which(str_detect(amniocytes.gestation.cell.line$variable, "^GW25"))] <- "Stage I (GW15-GW27)"
amniocytes.gestation.cell.line$stage[which(str_detect(amniocytes.gestation.cell.line$variable, "^GW26"))] <- "Stage I (GW15-GW27)"
amniocytes.gestation.cell.line$stage[which(str_detect(amniocytes.gestation.cell.line$variable, "^GW27"))] <- "Stage I (GW15-GW27)"
amniocytes.gestation.cell.line$stage[which(str_detect(amniocytes.gestation.cell.line$variable, "^GW28"))] <- "Stage II (GW28-GW30)"
amniocytes.gestation.cell.line$stage[which(str_detect(amniocytes.gestation.cell.line$variable, "^GW29"))] <- "Stage II (GW28-GW30)"
amniocytes.gestation.cell.line$stage[which(str_detect(amniocytes.gestation.cell.line$variable, "^GW30"))] <- "Stage II (GW28-GW30)"
p <- ggplot(data = amniocytes.gestation.cell.line, 
            mapping = aes(x = stage, y = value, fill = stage)) +
  geom_boxplot(width = 0.15, alpha = 0.9, outlier.shape = NA, color = "black", lwd = 0.4) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("#A63F94", "#F0C1DA")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0.05)) +
  stat_compare_means(
    method = "t.test",
    label = "p.format",
    label.x = 1.5,
    size = 4,
    vjust = 1.5
  ) +
  facet_grid(
    rows = vars(cell.line),
    scales = "free_y"
  ) +
  labs(x = "Stage", y = "Ratio") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    panel.grid = element_blank(),                # 移除网格线
    strip.background = element_rect(fill = "#f5f5f5", color = NA),
    strip.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )
ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/up.load/1.normal/5.0.output/5.0.amniocytes.cell.line.pvalue.20251005.pdf",width = 5,height = 8)

############################################################################################################################################################
###########################################celltype.ratio###################################################################################################
############################################################################################################################################################
table(metadata$singleR_labels)
amniocytes.celltype <- as.data.frame(table(metadata$singleR_labels))
colnames(amniocytes.celltype) <- c("Cell.type", "Freq")
amniocytes.celltype <- amniocytes.celltype %>%
  mutate(percent = Freq/sum(Freq))


amniocytes.celltype <- amniocytes.celltype[order(amniocytes.celltype$percent, decreasing = T),]
amniocytes.celltype$Cell.type <- factor(amniocytes.celltype$Cell.type, levels = amniocytes.celltype$Cell.type)
p <- ggplot(data=amniocytes.celltype)+
  geom_bar(stat = 'identity', mapping=aes(y = Cell.type, x = percent), fill = "#704D9E", alpha=0.8)+
  scale_x_continuous(expand = c(0,0), limits=c(0,0.30))+
  theme(axis.text.x = element_text(colour = "black", angle = 30, hjust = 1), 
        axis.text.y = element_text(colour = "black", hjust = 1),
        legend.position = "none", panel.background = element_rect(fill = 'white'), 
        axis.line.x = element_line(colour = "black",  size=0.5, lineend = "butt"), 
        axis.line.y = element_line(colour = "black", size=0.5)) +
  geom_text(
    mapping = aes(y = reorder(Cell.type, -percent), x = percent, label = paste0(round(percent*100, 2), "%")),
    stat = "identity",
    vjust = 0.5, hjust = -0.1)
ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/up.load/1.normal/5.0.output/5.0.amniocytes.celltype.ratio.20251005.pdf",width = 10,height = 8)

#######Dimplot by indivicual organ
col.Kidney <- c("#8CCAFD", "#009FF6")
col.Stomach <-c("#369891", "#86CB66","#91AD5A", "#007A0B")
col.Lung <- c( "#FF662E", "#FFC763")
col.Intestine <- c("#FC8DCA")
col.Placenta <- c( "#704D9E")
amniocytes.normal@meta.data$singleR_labels <- factor(amniocytes.normal@meta.data$singleR_labels, levels = c("Stomach-Myeloid.cells",                         
                                                                                                            "Stomach-Erythroblasts",                         
                                                                                                            "Stomach-Mesothelial.cells",                     
                                                                                                            "Stomach-Ciliated.epithelial.cells", 
                                                                                                            "Kidney-Megakaryocytes", 
                                                                                                            "Kidney-Ureteric.bud.cells", 
                                                                                                            "Lung-Squamous.epithelial.cells",                
                                                                                                            "Lung-Bronchiolar.and.alveolar.epithelial.cells",
                                                                                                            "Placenta-Megakaryocytes",                       
                                                                                                            "Intestine-Intestinal.epithelial.cells"))

color.cell.type <- c(col.Stomach, "#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE")
p <- DimPlot(object = amniocytes.normal, reduction = "umap", group.by = "singleR_labels", label = F, repel = TRUE, cols = color.cell.type)
ggsave (p,file= "/share/home/yangjingyi/project/1.BD/up.load/1.normal/5.0.output/5.0.dimplot.stomach.20251005.pdf",width = 8,height = 8)

###
amniocytes.normal@meta.data$singleR_labels <- factor(amniocytes.normal@meta.data$singleR_labels, levels = c("Kidney-Megakaryocytes", 
                                                                                                            "Kidney-Ureteric.bud.cells", 
                                                                                                            "Stomach-Myeloid.cells",                         
                                                                                                            "Stomach-Erythroblasts",                         
                                                                                                            "Stomach-Mesothelial.cells",                     
                                                                                                            "Stomach-Ciliated.epithelial.cells", 
                                                                                                            "Lung-Squamous.epithelial.cells",                
                                                                                                            "Lung-Bronchiolar.and.alveolar.epithelial.cells",
                                                                                                            "Placenta-Megakaryocytes",                       
                                                                                                            "Intestine-Intestinal.epithelial.cells"))

color.cell.type <- c(col.Kidney, "#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE")
p <- DimPlot(object = amniocytes.normal, reduction = "umap", group.by = "singleR_labels", label = F, repel = TRUE, cols = color.cell.type)
ggsave (p,file= "/share/home/yangjingyi/project/1.BD/up.load/1.normal/5.0.output/5.0.dimplot.kidney.20251005.pdf",width = 8,height = 8)

###
amniocytes.normal@meta.data$singleR_labels <- factor(amniocytes.normal@meta.data$singleR_labels, levels = c("Lung-Squamous.epithelial.cells",                
                                                                                                            "Lung-Bronchiolar.and.alveolar.epithelial.cells",
                                                                                                            "Kidney-Megakaryocytes", 
                                                                                                            "Kidney-Ureteric.bud.cells", 
                                                                                                            "Stomach-Myeloid.cells",                         
                                                                                                            "Stomach-Erythroblasts",                         
                                                                                                            "Stomach-Mesothelial.cells",                     
                                                                                                            "Stomach-Ciliated.epithelial.cells", 
                                                                                                            "Placenta-Megakaryocytes",                       
                                                                                                            "Intestine-Intestinal.epithelial.cells"))

color.cell.type <- c(col.Lung, "#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE")
p <- DimPlot(object = amniocytes.normal, reduction = "umap", group.by = "singleR_labels", label = F, repel = TRUE, cols = color.cell.type)
ggsave (p,file= "/share/home/yangjingyi/project/1.BD/up.load/1.normal/5.0.output/5.0.dimplot.lung.20251005.pdf",width = 8,height = 8)

###
amniocytes.normal@meta.data$singleR_labels <- factor(amniocytes.normal@meta.data$singleR_labels, levels = c("Placenta-Megakaryocytes",                       
                                                                                                            "Intestine-Intestinal.epithelial.cells",
                                                                                                            "Lung-Squamous.epithelial.cells",                
                                                                                                            "Lung-Bronchiolar.and.alveolar.epithelial.cells",
                                                                                                            "Kidney-Megakaryocytes", 
                                                                                                            "Kidney-Ureteric.bud.cells", 
                                                                                                            "Stomach-Myeloid.cells",                         
                                                                                                            "Stomach-Erythroblasts",                         
                                                                                                            "Stomach-Mesothelial.cells",                     
                                                                                                            "Stomach-Ciliated.epithelial.cells"))

color.cell.type <- c(col.Placenta,"#AEAEAE", "#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE")
p <- DimPlot(object = amniocytes.normal, reduction = "umap", group.by = "singleR_labels", label = F, repel = TRUE, cols = color.cell.type)
ggsave (p,file= "/share/home/yangjingyi/project/1.BD/up.load/1.normal/5.0.output/5.0.dimplot.plancenta.20251005.pdf",width = 8,height = 8)

###
amniocytes.normal@meta.data$singleR_labels <- factor(amniocytes.normal@meta.data$singleR_labels, levels = c("Intestine-Intestinal.epithelial.cells",
                                                                                                            "Lung-Squamous.epithelial.cells",                
                                                                                                            "Lung-Bronchiolar.and.alveolar.epithelial.cells",
                                                                                                            "Kidney-Megakaryocytes", 
                                                                                                            "Kidney-Ureteric.bud.cells", 
                                                                                                            "Stomach-Myeloid.cells",                         
                                                                                                            "Stomach-Erythroblasts",                         
                                                                                                            "Stomach-Mesothelial.cells",                     
                                                                                                            "Stomach-Ciliated.epithelial.cells","Placenta-Megakaryocytes"))

color.cell.type <- c(col.Intestine,"#AEAEAE", "#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE","#AEAEAE")
p <- DimPlot(object = amniocytes.normal, reduction = "umap", group.by = "singleR_labels", label = F, repel = TRUE, cols = color.cell.type)
ggsave (p,file= "/share/home/yangjingyi/project/1.BD/up.load/1.normal/5.0.output/5.0.dimplot.Intestine.20251005.pdf",width = 8,height = 8)


################################################################################################
##################pie.################################################################
################################################################################################
#############organ
metadata$organ <- NA
metadata$organ[which(str_detect(metadata$singleR_labels, "^Stomach"))] <- "Stomach"
metadata$organ[which(str_detect(metadata$singleR_labels, "^Kidney"))] <- "Kidney"
metadata$organ[which(str_detect(metadata$singleR_labels, "^Intestine"))] <- "Intestine"
metadata$organ[which(str_detect(metadata$singleR_labels, "^Lung"))] <- "Lung"
metadata$organ[which(str_detect(metadata$singleR_labels, "^Placenta"))] <- "Placenta"
#
mycolors <- c("#C8D589","#FAF3DD","#F9C382",
              "#C25CA8","#704D9E")

amniocytes.metadata.organ <- table(metadata$organ)
amniocytes.metadata.organ <- as.data.frame(amniocytes.metadata.organ)
colnames(amniocytes.metadata.organ) <- c("Organ", "Count")
amniocytes.metadata.organ <- amniocytes.metadata.organ %>%
  mutate(Ratio_Percent = scales::percent(Count / sum(Count)))


blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
amniocytes.metadata.organ <- amniocytes.metadata.organ[order(amniocytes.metadata.organ$Count, decreasing = T),]
amniocytes.metadata.organ$Organ <- factor(amniocytes.metadata.organ$Organ, levels = amniocytes.metadata.organ$Organ)
p <- ggplot(data=amniocytes.metadata.organ, mapping=aes(x="", y=Count, fill=Organ))+
  geom_bar(stat = "identity", width=0.5,position = 'stack')+
  scale_fill_manual(values = rev(mycolors))+
  coord_polar("y", start = 0)+
  blank_theme +
  theme(axis.text.x=element_blank()) 
ggsave(filename="5.0.organ.ratio.pie.20251005.pdf",p,width = 8, height = 8)
#################germ.layer
metadata$germ.layer <- NA
metadata$germ.layer[which(str_detect(metadata$singleR_labels, "^Intestine"))] <- "Endoderm"
metadata$germ.layer[which(str_detect(metadata$singleR_labels, "^Kidney"))] <- "Mesoderm"
metadata$germ.layer[which(str_detect(metadata$singleR_labels, "^Lung"))] <- "Endoderm"
metadata$germ.layer[which(str_detect(metadata$singleR_labels, "^Placenta"))] <- "Placenta"
metadata$germ.layer[which(str_detect(metadata$singleR_labels, "^Stomach"))] <- "Endoderm"

amniocytes.metadata.germ.layer <- table(metadata$germ.layer)
amniocytes.metadata.germ.layer <- as.data.frame(amniocytes.metadata.germ.layer)
colnames(amniocytes.metadata.germ.layer) <- c("Germ.layer", "Count")

amniocytes.metadata.germ.layer <- amniocytes.metadata.germ.layer %>%
  mutate(Ratio_Percent = scales::percent(Count / sum(Count)))



amniocytes.metadata.germ.layer <- amniocytes.metadata.germ.layer[order(amniocytes.metadata.germ.layer$Count, decreasing = T),]
amniocytes.metadata.germ.layer$Germ.layer <- factor(amniocytes.metadata.germ.layer$Germ.layer, 
                                                    levels = amniocytes.metadata.germ.layer$Germ.layer)

col.germ.layer <- c("#704D9E", "#C25CA8","#C8D589")
p <- ggplot(data=amniocytes.metadata.germ.layer, mapping=aes(x="", y=Count, fill=Germ.layer))+
  geom_bar(stat = "identity", width=0.5,position = 'stack')+
  scale_fill_manual(values = alpha(col.germ.layer, .8))+
  coord_polar("y", start = 0)+
  blank_theme +
  theme(axis.text.x=element_blank()) 
ggsave(filename="5.0.germ.layer.ratio.pie.20251005.pdf",p,width = 8, height = 8)
