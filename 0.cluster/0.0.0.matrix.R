rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/1.cluster/0.0.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") #install


library(loomR)
lfile <- connect(filename = "/share/home/yangjingyi/project/1.BD/up.load/1.cluster/0.0.output/GSE156793_S3_gene_count.loom", mode = "r+")

for (i in 0:405){
  m = 10000 * i + 1
  n = 10000 * (i + 1)
  
  matrix2=lfile[["matrix"]][m:n,]
  matrix2=t(matrix2)
  colnames(matrix2)=barcode[m:n]
  row.names(matrix2)=gene[]
  print(i)
  write.csv(x = matrix2, row.names = TRUE, file = paste0("/share/home/yangjingyi/project/fetal_cell_sciencepaper/1.0amniocytes_reference_fetal_cells/",i,"_reference.csv"))
}

##
i=406
m = 10000 * i + 1
matrix2=lfile[["matrix"]][m:length(barcode),]
matrix2=t(matrix2)
colnames(matrix2)=barcode[m:length(barcode)]
row.names(matrix2)=gene[]
print(i)
write.csv(x = matrix2, row.names = TRUE, file = paste0("/share/home/yangjingyi/project/fetal_cell_sciencepaper/1.0amniocytes_reference_fetal_cells/",i,"_reference.csv"))


###
file[["col_attrs"]]
cell_type = lfile$col.attrs$Main_cluster_name[]
length(cell_type)

barcode=lfile$col.attrs$obs_names[]
length(barcode)
Organ_cell_lineage = lfile$col.attrs$Organ_cell_lineage$read()
organ = lfile$col.attrs$Organ[]

reference_data <- data.frame(barcode, Organ_cell_lineage, organ)
reference_data[,4] <- rownames(reference_data)

write.csv(x = reference_data, row.names = TRUE, file = "reference_data.csv")

