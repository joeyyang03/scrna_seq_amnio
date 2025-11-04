rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/2.model/3.1.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 

library(dplyr)
library(readr)
library(Seurat)
library(pls)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(tidyr)
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/0.0.output/amnio.train.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/3.0.output/3.0.amniocytes.with.celltype.cor.top30.exp.0.1.df.20251007.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/0.0.output/amnio.test.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/3.0.output/3.0.amniocytes.test.with.celltype.cor.top30.exp.0.1.df.20251007.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/5.0.output/amniocytes.cohort1.FGR.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/3.0.output/3.0.amniocytes.FGR.cohort1.with.celltype.cor.top30.exp.0.1.df.20251007.Robj")
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
p <- list()
p.val <- list()
for (i in 1:length(organ)){
  m <- organ[[i]]
  train.data <-  df[, grepl(m, names(df)) | grepl("orig\\.ident", names(df))]
  train.data<- merge(train.data,gestation.day.train,by="orig.ident")
  rownames(train.data) <- train.data$orig.ident
  train.data <- train.data[,-1]
  train.ges.day <- train.data[,c("gestation.day"), drop = FALSE]
  
  x <- model.matrix(gestation.day~.,train.data)[,-1]
  y <- train.data$gestation.day
  pls_model <- plsr(y ~ x, data = train.data, ncomp = 3,  validation = "CV", segments = 10) 
  
  
  ####
  pred_values <- predict(pls_model, x)[, 1, 3]  
  results <- data.frame(
    orig.ident = rownames(train.data),
    Chronological.age = y,
    Predicted.age = pred_values,
    group = "Normal.train"
  )
  
  ##############
  test.data <-  df.test[, grepl(m, names(df.test)) | grepl("orig\\.ident", names(df.test))]
  test.data <- merge(test.data,gestation.day.test,by="orig.ident")
  rownames(test.data) <- test.data$orig.ident
  test.data <- test.data[,-1]
  
  
  x.test <- model.matrix(gestation.day~.,test.data)[,-1]
  test.predictions <- predict(pls_model, newdata = x.test)
  
  
  test.predictions <- as.data.frame(test.predictions[,,1] )
  colnames(test.predictions) <- c("Predicted.age")
  test.predictions$orig.ident <- rownames(test.predictions)
  test.predictions <- merge(test.predictions,gestation.day.test,by="orig.ident")
  colnames(test.predictions) <- c("orig.ident", "Predicted.age", "Chronological.age")
  test.predictions$group <- "Normal.test"
  matrix <- rbind(results,test.predictions)
  ########
  p[[i]] <- ggscatter(results, x = "Chronological.age", y = "Predicted.age",          
                  add = "reg.line",           
                  conf.int = TRUE,    
                  conf.int.level = 0.99,
                  add.params = list(color = "#254689",                             
                                    fill = "#a0c9e5"),        
                  cor.coef = TRUE,#添加相关系数
                  cor.coeff.args = list(method = "pearson", label.x = 120, label.sep = "\n"))+
    theme(panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank()   
    )+ggtitle(paste0("Discovery cohort ",  m,"(Pls model)"))
  ###
  p.val[[i]] <- ggscatter(test.predictions, x = "Chronological.age", y = "Predicted.age",          
                  add = "reg.line",           
                  conf.int = TRUE,    
                  conf.int.level = 0.99,
                  add.params = list(color = "#B4355E",                             
                                    fill = "#f3dee0"),          
                  cor.coef = TRUE,
                  cor.coeff.args = list(method = "pearson", label.x = 120, label.sep = "\n"))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()   
    )+ggtitle(paste0("Validation cohort ",m,"(Pls model)"))
  
  #####################
  
  FGR.cohort1.data <-  df.FGR.cohort1[, grepl(m, names(df.FGR.cohort1)) | grepl("orig\\.ident", names(df.FGR.cohort1))]
  FGR.cohort1.data <- merge(FGR.cohort1.data,gestation.day.FGR.cor1,by="orig.ident")
  rownames(FGR.cohort1.data) <- FGR.cohort1.data$orig.ident
  FGR.cohort1.data <- FGR.cohort1.data[,-1]
  
  x.FGR.cohort1 <- model.matrix(gestation.day~.,FGR.cohort1.data)[,-1]
  FGR.cohort1.predictions <- predict(pls_model, newdata = x.FGR.cohort1)
  FGR.cohort1.predictions <- as.data.frame(FGR.cohort1.predictions[, , 1] )
  colnames(FGR.cohort1.predictions) <- c("Predicted.age")
  FGR.cohort1.predictions$orig.ident <- rownames(FGR.cohort1.predictions)
  FGR.cohort1.predictions <- merge(FGR.cohort1.predictions,gestation.day.FGR.cor1,by="orig.ident")
  colnames(FGR.cohort1.predictions) <- c("orig.ident", "Predicted.age", "Chronological.age")
  FGR.cohort1.predictions$group <- "FGR.cohort1"
  
  
  ##########################################################plot
  matrix <- rbind(matrix,FGR.cohort1.predictions)
  matrix$delta <- matrix$Chronological.age-matrix$Predicted.age
  
  matrix$organ <- m
  organ.results <- rbind(organ.results,matrix)
  print(m)
}
plot <- cowplot::plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],ncol = 5)
ggsave (plot,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/3.1.output/3.1.ggscatter.pls.model.cor.top30.exp.0.1.train.data.by.pls.20251010.pdf",width = 16,height = 4)
plot <- cowplot::plot_grid(p.val[[1]],p.val[[2]],p.val[[3]],p.val[[4]],p.val[[5]],ncol = 5)
ggsave (plot,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/3.1.output/3.1.ggscatter.pls.model.cor.top30.exp.0.1.test.data.by.pls.20251010.pdf",width = 16,height = 4)



res.sel <- organ.results[organ.results$group == "FGR.cohort1",]
unique(res.sel$orig.ident)

res.sel.plot <- res.sel %>%
  pivot_wider(
    id_cols = c(orig.ident, Chronological.age),  
    names_from = organ,          
    values_from = delta,         
  ) %>%
  as.data.frame()  

cols_to_scale <- 3:7
df_scaled <- res.sel.plot
df_scaled[, cols_to_scale] <- scale(res.sel.plot[, cols_to_scale], center = F)
df_long <- pivot_longer(df_scaled,
                        cols = Intestine:Stomach,
                        names_to = "Organ",
                        values_to = "Value")

p <- ggplot(df_long, aes(x = Chronological.age, y = Value, color = Organ)) +
  geom_point(alpha = 0.5, size = 2) +   
  coord_cartesian(ylim = c(-2, 2) )+
  geom_smooth(method = "loess", se = FALSE, size = 1.2, span = 0.75) + 
  labs(
    x = "Chronological age (days)",
    y = "Scaled expression",
    color = "Organ"
  ) +
  scale_color_manual(values = c(
    "Intestine" = "#B366C7",  
    "Kidney"    = "#F8DF79",  
    "Lung"      = "#E9212C",  
    "Placenta"  = "#414986",  
    "Stomach"   = "#FF8C42"   
  )) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),                
    panel.background = element_rect(fill = "white", color = NA), 
    axis.line = element_line(color = "black", linewidth = 0.6),  
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11)
  )



ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/3.1.output/3.1.ggplot.organ.scaled.stageII.stageIII.by.pls.20251007.pdf",width = 8,height = 8)
