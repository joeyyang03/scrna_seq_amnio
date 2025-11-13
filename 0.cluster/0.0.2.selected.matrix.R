rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/1.cluster/0.0.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") #install

library(dplyr)

#####
reference_data = read.csv(file = "/share/home/yangjingyi/project/1.BD/up.load/1.cluster/0.0.output/reference_172_celltype_based_1000.csv", header=T, sep=",",stringsAsFactors = F,row.names = 1)
reference_data$barcode
reference_data$barcode <- gsub('[-]', '.', reference_data$barcode)

###########
reference_matrix_selected <- read.csv("/share/home/yangjingyi/project/fetal_cell_sciencepaper/1.0amniocytes_reference_fetal_cells/0_reference.csv",header = T,stringsAsFactors = F)
reference_matrix_selected[1:10,1:10]
table(duplicated(reference_matrix_selected$X))
reference_matrix_selected=distinct(reference_matrix_selected, reference_matrix_selected$X, .keep_all = T)
rownames(reference_matrix_selected) <- reference_matrix_selected$X
table(colnames(reference_matrix_selected) %in% reference_data$barcode)
reference_matrix_selected <- reference_matrix_selected[,colnames(reference_matrix_selected) %in% reference_data$barcode]
reference_matrix_selected$gene <- rownames(reference_matrix_selected)

for (i in 1:406){
  df <- read.csv(paste0("/share/home/yangjingyi/project/fetal_cell_sciencepaper/1.0amniocytes_reference_fetal_cells/",i,"_reference.csv"),header = T,stringsAsFactors = F)
  df[1:10,1:10]
  table(duplicated(df$X))
  df=distinct(df, df$X, .keep_all = T)
  rownames(df) <- df$X
  table(colnames(df) %in% reference_data$barcode)
  df.selected <- df[,colnames(df) %in% reference_data$barcode]
  df.selected$gene <- rownames(df.selected)
  reference_matrix_selected <- merge(reference_matrix_selected, df.selected,by="gene")
  
  print(i)
}


write.csv(x = reference_matrix_selected, row.names = TRUE, file = "/share/home/yangjingyi/project/1.BD/up.load/1.cluster/0.0.output/reference_selected_matrix_20250326.csv", sep=",")
