rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/2.model/2.1.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 

library(dplyr)
library(readr)
library(Seurat)
library(caret)
library(randomForest)
library(ggplot2)
library(ggpubr)
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/0.0.output/amnio.train.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/2.0.output/2.0.amniocytes.with.celltype.cor.0.8.exp.0.1.df.20251007.Robj")
##################################################
gestation.day.train <- amniocytes.sample@meta.data
gestation.day.train <- gestation.day.train[,c("gestation.day", "orig.ident")]
gestation.day.train <- gestation.day.train %>%
  group_by(orig.ident) %>%
  summarise(
    gestation.day = first(gestation.day)
  )

df <- merge(df,gestation.day.train,by="orig.ident")
##################################################
train.data <- df
train.data[is.na(train.data)] = 0
rownames(train.data) <- train.data$orig.ident
train.data <- train.data[,-1]
colnames(train.data) <- gsub("-",".",colnames(train.data))
##########################################################################################################################
set.seed(123)
ctrl <- rfeControl(functions = rfFuncs, method = "cv", number = 5)
ncol(train.data)
#[1] 974
rfe_model <- rfe(
  x = train.data[, -ncol(train.data)],
  y = train.data$gestation.day,
  sizes = 1:(ncol(train.data)-1),
  rfeControl = ctrl,
  metric = "RMSE" 
)
save(rfe_model, file = "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.1.output/2.1.rfe_model.20251010.Robj")

print(rfe_model)
rfe_res <- rfe_model$results
p <- ggplot(rfe_res, aes(x = Variables, y = Rsquared)) +
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
  geom_vline(xintercept = 10, linetype = "dashed", color = "steelblue") +
  geom_vline(xintercept = 100, linetype = "dashed", color = "brown") +
  labs(x = "Number of Features", y = "Rsquared")
ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.1.output/2.1.features.number.selection.2025.pdf",width = 40,height = 8)

p <- ggplot(rfe_res[1:100,], aes(x = Variables, y = Rsquared)) +
  geom_line(size = 0.5, color = "#af93c4", alpha = 0.5) +
  geom_point(size = 1, color = "#A63F94", alpha = 0.8) +
  theme(
    panel.grid = element_blank(),     
    panel.border = element_blank(),    
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.background = element_rect(fill = "white"),
    axis.title = element_text(face = "bold")
  ) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "steelblue", size =1) +
  labs(x = "Number of Features", y = "Rsquared")
ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.1.output/2.1.features.number.selection.top100.2025.pdf",width = 8,height = 8)
