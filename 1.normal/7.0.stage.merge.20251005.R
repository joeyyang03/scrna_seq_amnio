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

#############################################################################
#############################################################################
#############################################################################
expr <- amniocytes.normal[["RNA"]]@data  

table(metadata$gestation.week,metadata$orig.ident)
group_info <- paste(metadata$gestation.week, metadata$singleR_labels, sep = "__")

avg_exp <- sapply(unique(group_info), function(g) {
  cells <- which(group_info == g)
  if (length(cells) == 1) {
    expr[, cells]
  } else {
    Matrix::rowMeans(expr[, cells, drop = FALSE])
  }
})

avg_exp <- as.data.frame(avg_exp)
colnames(avg_exp) <- unique(group_info)
avg_exp[1:3,1:3]

avg_exp$gene <- rownames(expr)  



avg_exp_wide <- avg_exp %>%
  tidyr::pivot_longer(-gene, names_to = "sample_celltype", values_to = "expr") %>%
  tidyr::separate(sample_celltype, into = c("sample", "celltype"), sep = "__", extra = "merge") %>%
  tidyr::unite("celltype_gene", celltype, gene, sep = "__", remove = FALSE) %>%  # 保留 gene 防止丢失
  tidyr::pivot_wider(names_from = sample, values_from = expr) %>%
  tibble::column_to_rownames("celltype_gene")

avg_exp_wide <- avg_exp_wide %>% 
  dplyr::select(-gene, -celltype)

table(is.na(avg_exp_wide))


avg_exp_wide[is.na(avg_exp_wide)] <- 0
dim(avg_exp_wide)


filtered_mat <- avg_exp_wide[rowSums(avg_exp_wide > 1) >= 2, ]
dim(filtered_mat)


cv_values <- apply(filtered_mat, 1, function(x) {
  if (mean(x) == 0) {
    return(NA)  
  } else {
    return(sd(x) / mean(x))
  }
})


filtered_cv <- filtered_mat[cv_values > 0.5, ]
dim(filtered_cv)

cor(filtered_cv)
pheatmap(cor(filtered_cv),
         color = colorRampPalette(c("#F5A551", "white", "#A277AF"))(50), filename = "7.0.cor.gestation.week.matrix.20251005.pdf")




save.image(file = "/share/home/yangjingyi/project/1.BD/up.load/1.normal/7.0.output/7.0.cor.gestation.week.matrix.20251005.RData")


