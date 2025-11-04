rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/1.cluster/0.0.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 
library(Seurat)
library(rhdf5)
library("DropletUtils")

folder_path <- "/share/home/yangjingyi/project/1.BD/data/rawdata.20251001/"

dir <- list.files(folder_path, full.names = FALSE, recursive = FALSE)
samples <- grep(".h5$", dir, value = TRUE)

print(samples)

for (i in 1:length(samples)){
  sample <- samples[[i]]
  dir <- paste0("/share/home/yangjingyi/project/1.BD/data/rawdata.20251001/", sample)
  amniocyte.info = h5read(dir,"/")
  amniocyte.info1 = downsampleReads(dir,prop = 1)
  
  
  geneid = amniocyte.info$features$id
  genename = amniocyte.info$features$name
  cellname = paste(amniocyte.info$barcodes[amniocyte.info$barcode_info$pass_filter[1,]+1],"-","1",sep = "")
  filter.matrix = amniocyte.info1[,cellname]
  
  prefixes = gsub(".h5","",sample)
  write10xCounts(paste0("/share/home/yangjingyi/project/1.BD/data/rawdata.20251001/",prefixes,".matrix"),filter.matrix,gene.id=geneid,gene.symbol=genename,version="3")
  
}