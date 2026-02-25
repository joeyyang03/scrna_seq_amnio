rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/2.model/3.1.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 

library(caret)
library(randomForest)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(tidyr)
library(dplyr)
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
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(
    gestation.day = dplyr::first(gestation.day)
  )

df[is.na(df)] = 0
#################################################
gestation.day.test <- amniocytes.test@meta.data
gestation.day.test <- gestation.day.test[,c("gestation.day", "orig.ident")]
gestation.day.test <- gestation.day.test %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(
    gestation.day = dplyr::first(gestation.day)
  )
df.test[is.na(df.test)] = 0
##################################################
gestation.day.FGR.cor1 <- amniocytes.FGR.cohort1@meta.data
gestation.day.FGR.cor1 <- gestation.day.FGR.cor1[,c("gestation.day", "orig.ident")]
gestation.day.FGR.cor1 <- gestation.day.FGR.cor1 %>%
  dplyr::group_by(orig.ident) %>%
  dplyr::summarise(
    gestation.day = dplyr::first(gestation.day)
  )
df.FGR.cohort1[is.na(df.FGR.cohort1)] = 0


######################################################
organ.results <- data.frame(matrix(nrow=0,ncol=6))
colnames(organ.results) <- c("Chronological.age", "Predicted.age", "group" , "delta", "organ", "orig.ident")

organ <- c("Intestine",    "Kidney",      "Lung",  "Placenta",   "Stomach")

rfe.variant.num <- data.frame(organ = c("Intestine",    "Kidney",      "Lung",  "Placenta",   "Stomach"),
                              num = c(13,25,11,10,15))
p.variant.num <- list()
p.variant.cor <- list()
p <- list()
p.val <- list()
for (i in 1:length(organ)){
  m <- organ[[i]]
  train.data <-  df[, grepl(m, names(df)) | grepl("orig\\.ident", names(df))]
  train.data<- merge(train.data,gestation.day.train,by="orig.ident")
  rownames(train.data) <- train.data$orig.ident
  train.data <- train.data[,-1]
  colnames(train.data) <- gsub("-",".",colnames(train.data))
  
  N <- rfe.variant.num[rfe.variant.num$organ == m,]$num#25

  set.seed(123)

  train_control <- trainControl(
    method = "cv",
    number = 5
  )
  
  rf.model <- train(
    x = train.data[, -ncol(train.data)],
    y = train.data$gestation.day,
    method = "rf",
    trControl = train_control,
    importance = TRUE
  )
save(rf.model, file = paste0("/share/home/yangjingyi/project/1.BD/up.load/2.model/3.1.output/3.1.rf.model.",m,".20251010.Robj"))
  
  var_imp <- varImp(rf.model, scale = TRUE)
  print(var_imp)
  plot(var_imp)

  imp_df <- var_imp$importance
  
  imp_df$Variable <- rownames(imp_df)
  
  imp_df_sorted <- imp_df[order(imp_df$Overall, decreasing = T), ]
  imp_df_sorted.plot <- imp_df_sorted[1:N,]
  imp_df_sorted.plot <- imp_df_sorted.plot[order(imp_df_sorted.plot$Overall, decreasing = F), ]
  imp_df_sorted.plot$Variable <- factor(imp_df_sorted.plot$Variable, levels = imp_df_sorted.plot$Variable)
  p.variant.cor[[i]] <- ggplot(data=imp_df_sorted.plot)+
    geom_bar(stat = 'identity', mapping=aes(y = Variable, x = Overall), fill = "#af93c4")+
    scale_x_continuous(expand = c(0,0), limits=c(0,100))+
    theme(axis.text.x = element_text(colour = "black", angle = 30, hjust = 1), 
          axis.text.y = element_text(colour = "black", hjust = 1),
          legend.position = "none", panel.background = element_rect(fill = 'white'), 
          axis.line.x = element_line(colour = "black",  size=0.5, lineend = "butt"), 
          axis.line.y = element_line(colour = "black", size=0.5)) 
 ggsave (p.variant.cor[[i]] ,file= paste0("/share/home/yangjingyi/project/1.BD/up.load/2.model/3.1.output/3.1.importance.N.variable.rf.select.",m,".20251007.pdf"),width = 8,height = 8)
  ####
  top_vars <- imp_df_sorted.plot$Variable
  train_control <- trainControl(method = "cv", number = 5, savePredictions = TRUE)
  
  pls_model <- train(
    as.formula(paste("gestation.day ~", paste(top_vars, collapse = " + "))),
    data = train.data,
    method = "pls",
    trControl = train_control,
    tuneLength = 10,
    preProcess = c("center", "scale")
  )
  save(pls_model, file = paste0("/share/home/yangjingyi/project/1.BD/up.load/2.model/3.1.output/3.1.pls_model.rf.select.",m,".20251010.Robj"))
  #####plot
  best_ncomp <- pls_model$bestTune$ncomp
  pred_train <- pls_model$pred
  pred_train_best <- pred_train %>%
    dplyr::filter(ncomp == best_ncomp)
  
  train_pred_df <- train.data %>%
    dplyr::mutate(
      Predicted_Gestation.day = pred_train_best$pred[match(
        1:nrow(train.data),
        pred_train_best$rowIndex
      )]
    )
  #######
  p[[i]] <- ggscatter(train_pred_df, x = "gestation.day", y = "Predicted_Gestation.day",          
                      add = "reg.line",           
                      conf.int = TRUE,    
                      conf.int.level = 0.99,
                      add.params = list(color = "#254689",                             
                                        fill = "#a0c9e5"),         
                      cor.coef = TRUE,
                      cor.coeff.args = list(method = "pearson", label.x = 120, label.sep = "\n"))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank() 
    )+coord_cartesian(xlim = c(100, 230), ylim = c(100, 230)) +
    ggtitle(paste0("Discovery cohort ",  m,"(Pls model)"))

  ##############test
  test.data <-  df.test[, grepl(m, names(df.test)) | grepl("orig\\.ident", names(df.test))]
  test.data <- merge(test.data,gestation.day.test,by="orig.ident")
  rownames(test.data) <- test.data$orig.ident
  test.data <- test.data[,-1]
  colnames(test.data) <- gsub("-",".",colnames(test.data))
  
  pred_test <- predict(
    pls_model,
    newdata = test.data
  )
  test.data <- test.data %>%
    dplyr::mutate(
      Predicted_Gestation.day = as.numeric(pred_test)
    )
  cor(test.data$gestation.day, test.data$Predicted_Gestation.day)

  #######
  p.val[[i]] <- ggscatter(test.data, x = "gestation.day", y = "Predicted_Gestation.day",          
                          add = "reg.line",           
                          conf.int = TRUE,    
                          conf.int.level = 0.99,
                          add.params = list(color = "#B4355E",                             
                                            fill = "#f3dee0"),       
                          cor.coef = TRUE,#添加相关系数
                          cor.coeff.args = list(method = "pearson", label.x = 140, label.sep = "\n"))+
    theme(panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank()  
    )+
    coord_cartesian(xlim = c(110, 230), ylim = c(110, 230)) +
    ggtitle(paste0("Validation cohort ",m,"(Pls model)"))

  #####################
  
  FGR.cohort1.data <-  df.FGR.cohort1[, grepl(m, names(df.FGR.cohort1)) | grepl("orig\\.ident", names(df.FGR.cohort1))]
  FGR.cohort1.data <- merge(FGR.cohort1.data,gestation.day.FGR.cor1,by="orig.ident")
  rownames(FGR.cohort1.data) <- FGR.cohort1.data$orig.ident
  FGR.cohort1.data <- FGR.cohort1.data[,-1]
  colnames(FGR.cohort1.data) <- gsub("-",".",colnames(FGR.cohort1.data))
  
  pred_FGR.cohort1 <- predict(
    pls_model,
    newdata = FGR.cohort1.data
  )
  FGR.cohort1.data <- FGR.cohort1.data %>%
    dplyr::mutate(
      Predicted_Gestation.day = as.numeric(pred_FGR.cohort1)
    )
  ##########################################################plot
  train_pred_df$group <- "Train"
  test.data$group <- "Test"
  FGR.cohort1.data$group <- "FGR.cohort1"
  matrix <- rbind(train_pred_df[,c("gestation.day","Predicted_Gestation.day","group")],
                  test.data[,c("gestation.day","Predicted_Gestation.day","group")],
                  FGR.cohort1.data[,c("gestation.day","Predicted_Gestation.day","group")])
  matrix$delta <- matrix$gestation.day-matrix$Predicted_Gestation.day
  
  matrix$organ <- m
  matrix$orig.ident <- rownames(matrix)
  colnames(matrix) <- c("Chronological.age", "Predicted.age", "group" , "delta", "organ", "orig.ident")
  organ.results <- rbind(organ.results,matrix)
  print(m)
  save.image(file = paste0("/share/home/yangjingyi/project/1.BD/up.load/2.model/3.1.output/3.1.pls_model.",m,".20251010.RData"))
  
}


###################
plot <- cowplot::plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],ncol = 5)
ggsave (plot,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/3.1.output/3.1.ggscatter.pls.model.cor.0.8.exp.0.1.train.data.by.pls.20251010.pdf",width = 16,height = 4)
plot <- cowplot::plot_grid(p.val[[1]],p.val[[2]],p.val[[3]],p.val[[4]],p.val[[5]],ncol = 5)
ggsave (plot,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/3.1.output/3.1.ggscatter.pls.model.cor.0.8.exp.0.1.test.data.by.pls.20251010.pdf",width = 16,height = 4)
res.sel <- organ.results[organ.results$group == "FGR.cohort1",]
unique(res.sel$orig.ident)


res.sel.plot <- res.sel %>%
  pivot_wider(
    id_cols    =  c(orig.ident, Chronological.age),
    names_from = organ,
    values_from = delta
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



ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/3.1.output/3.1.ggplot.cor.0.8.organ.scaled.stageII.stageIII.by.pls.20251010.pdf",width = 8,height = 8)


