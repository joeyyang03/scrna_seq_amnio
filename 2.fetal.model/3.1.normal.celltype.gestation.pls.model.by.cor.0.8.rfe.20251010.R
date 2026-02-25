rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/2.model/3.1.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 

library(dplyr)
library(readr)
library(Seurat)
library(caret)
library(randomForest)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(tidyr)
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/0.0.output/amnio.train.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/3.0.output/3.0.amniocytes.with.celltype.cor.without.selection.exp.0.1.df.20251007.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/0.0.output/amnio.test.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/3.0.output/3.0.amniocytes.test.with.celltype.cor.without.selection.exp.0.1.df.20251007.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/5.0.output/amniocytes.cohort1.FGR.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/3.0.output/3.0.amniocytes.FGR.cohort1.with.celltype.cor.without.selection.exp.0.1.df.20251007.Robj")
##################################################
gestation.day.train <- amniocytes.sample@meta.data
gestation.day.train <- gestation.day.train[,c("gestation.day", "orig.ident")]
gestation.day.train <- gestation.day.train %>%
  group_by(orig.ident) %>%
  summarise(
    gestation.day = first(gestation.day)
  )

df[is.na(df)] = 0
#################################################
gestation.day.test <- amniocytes.test@meta.data
gestation.day.test <- gestation.day.test[,c("gestation.day", "orig.ident")]
gestation.day.test <- gestation.day.test %>%
  group_by(orig.ident) %>%
  summarise(
    gestation.day = first(gestation.day)
  )
df.test[is.na(df.test)] = 0
##################################################
gestation.day.FGR.cor1 <- amniocytes.FGR.cohort1@meta.data
gestation.day.FGR.cor1 <- gestation.day.FGR.cor1[,c("gestation.day", "orig.ident")]
gestation.day.FGR.cor1 <- gestation.day.FGR.cor1 %>%
  group_by(orig.ident) %>%
  summarise(
    gestation.day = first(gestation.day)
  )
df.FGR.cohort1[is.na(df.FGR.cohort1)] = 0


######################################################
organ.results <- data.frame(matrix(nrow=0,ncol=6))
colnames(organ.results) <- c("orig.ident",       "Predicted.age",     "Chronological.age", "group" , "delta", "organ")

organ <- c("Intestine",    "Kidney",      "Lung",  "Placenta",   "Stomach")

p.variant.num <- list()
for (i in 1:length(organ)){
  m <- organ[[i]]
  train.data <-  df[, grepl(m, names(df)) | grepl("orig\\.ident", names(df))]
  train.data<- merge(train.data,gestation.day.train,by="orig.ident")
  rownames(train.data) <- train.data$orig.ident
  train.data <- train.data[,-1]
  colnames(train.data) <- gsub("-",".",colnames(train.data))

  set.seed(123)
  ctrl <- rfeControl(functions = rfFuncs, method = "cv", number = 5)
  ncol(train.data)
  
  rfe_model <- rfe(
    x = train.data[, -ncol(train.data)],
    y = train.data$gestation.day,
    sizes = 1:(ncol(train.data)-1),
    rfeControl = ctrl,
    metric = "RMSE" 
  )
  print(rfe_model)
  rfe_res <- rfe_model$results
  p.variant.num[[i]] <- ggplot(rfe_res, aes(x = Variables, y = Rsquared)) +
    geom_line(size = 1, color = "#af93c4", alpha = 0.5) +
    geom_point(size = 1.5, color = "#A63F94", alpha = 0.8) +
    theme(
      panel.grid = element_blank(),     
      panel.border = element_blank(),  
      axis.line = element_line(color = "black"), 
      axis.ticks = element_line(color = "black"),
      panel.background = element_rect(fill = "white"),
      axis.title = element_text(face = "bold")
    ) +
    labs(x = "Number of Features", y = "Rsquared")
  ggsave (p.variant.num[[i]],file= paste0("/share/home/yangjingyi/project/1.BD/up.load/2.model/3.1.output/3.1.feature.",m,".20251010.pdf"),width = 45,height = 4)
  
  p.variant.num[[i]] <- ggplot(rfe_res[1:100,], aes(x = Variables, y = Rsquared)) +
    geom_line(size = 1, color = "#af93c4", alpha = 0.5) +
    geom_point(size = 1.5, color = "#A63F94", alpha = 0.8) +
    theme(
      panel.grid = element_blank(),  
      panel.border = element_blank(),   
      axis.line = element_line(color = "black"), 
      axis.ticks = element_line(color = "black"),
      panel.background = element_rect(fill = "white"),
      axis.title = element_text(face = "bold")
    ) +
    labs(x = "Number of Features", y = "Rsquared")
  ggsave (p.variant.num[[i]],file= paste0("/share/home/yangjingyi/project/1.BD/up.load/2.model/3.1.output/3.1.feature.",m,".top100.20251010.pdf"),width = 8,height = 8)
}