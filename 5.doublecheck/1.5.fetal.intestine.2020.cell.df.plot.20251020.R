rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/5.doublecheck/1.5.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") #install
library(hdf5r)
library(Matrix)
library(Seurat)
library(tidyverse)
library(dplyr)
library(tidyr)
amnio.matrix <- read.csv("/share/home/yangjingyi/project/1.BD/up.load/5.doublecheck/1.5.output/1.5.amnio.matrix.2020.cell.epi.20251025.csv", row.names = 1)
fetal.matrix <- read.csv("/share/home/yangjingyi/project/1.BD/up.load/5.doublecheck/1.5.output/1.5.fetal.matrix.2020.cell.epi.20251025.csv", row.names = 1)

summ.matrix <-  data.frame(matrix(ncol = 4, nrow = 0)) 
colnames(summ.matrix) <- c("gestation.week", "num.fetal", "num.amnio", "num.overlap")

ges.week <- c("GW17", "GW18", "GW19", "GW20", "GW21", "GW22", "GW24")
for (i in 1:length(ges.week)){
  ges <- ges.week[[i]]
  amnio.sel <- amnio.matrix[amnio.matrix$gestation.week %in% ges,]
  rownames(amnio.sel) <- amnio.sel$gestation.week
  amnio.sel <- amnio.sel %>% select(-gestation.week)
  amnio.sel <- amnio.sel[,amnio.sel > 0.1]
  
  fetal.sel <- fetal.matrix[fetal.matrix$gestation.week %in% ges,]
  rownames(fetal.sel) <- fetal.sel$gestation.week
  fetal.sel <- fetal.sel %>% select(-gestation.week)
  fetal.sel <- fetal.sel[,fetal.sel > 0.1]
  
  matrix.tmp <-  data.frame(matrix(ncol = 4, nrow = 1)) 
  colnames(matrix.tmp) <- c("gestation.week", "num.fetal", "num.amnio", "num.overlap")
  matrix.tmp$gestation.week <- ges
  matrix.tmp$num.overlap <- length(intersect(colnames(amnio.sel), colnames(fetal.sel)))
  matrix.tmp$num.fetal <- ncol(fetal.sel) - length(intersect(colnames(amnio.sel), colnames(fetal.sel)))
  matrix.tmp$num.amnio <- ncol(amnio.sel) - length(intersect(colnames(amnio.sel), colnames(fetal.sel)))

  summ.matrix <- rbind(summ.matrix,matrix.tmp)
}


df_long <- summ.matrix %>%
  pivot_longer(
    cols = c(num.fetal, num.amnio, num.overlap),
    names_to = "type",
    values_to = "count"
  )

df_long$type <- factor(df_long$type, levels = c("num.fetal", "num.overlap", "num.amnio"))


p <- ggplot(df_long, aes(x = gestation.week, y = count, fill = type)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = count), color = "black", size = 3,
            position = position_stack(vjust = 0.5)) + 
  scale_fill_manual(values = c("num.fetal" = "#f3dee0",
                               "num.amnio" = "#a0c9e5",
                               "num.overlap" = "#817cb9")) +
  labs(x = "Gestation Week", y = "Count", fill = "Type") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.7),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+ggtitle("Intestine")
ggsave(filename="1.5.gene.num.ges.20251026.pdf",p,width = 8, height = 8)
