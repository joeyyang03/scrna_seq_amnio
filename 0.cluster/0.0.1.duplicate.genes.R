rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/1.cluster/0.0.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") #install

library(dplyr)

for (i in 0:406){
    df <- read.csv(paste0("/share/home/yangjingyi/project/fetal_cell_sciencepaper/1.0amniocytes_reference_fetal_cells/",i,"_reference.csv"),header = T)
    head(df)
    table(duplicated(df$X))
    df=distinct(df, df$X, .keep_all = T)
    write.csv(x = df, row.names = TRUE, file = paste0("/share/home/yangjingyi/project/1.BD/up.load/1.cluster/0.0.output/clear_",i,"_reference.csv"))
    print(i)
    }