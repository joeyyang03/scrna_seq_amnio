rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/2.model/2.1.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") #install

library(dplyr)
library(caret)
library(randomForest)
library(ggplot2)
library(ggpubr)
library(pROC)
library(stringr)
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/0.0.output/amnio.train.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/2.0.output/2.0.amniocytes.with.celltype.cor.0.8.exp.0.1.df.20251007.Robj")
##################################################
gestation.day.train <- amniocytes.sample@meta.data
gestation.day.train <- gestation.day.train[,c("gestation.day", "orig.ident")]
gestation.day.train <- gestation.day.train %>%
  group_by(orig.ident) %>%
  summarise(
    gestation.day = dplyr::first(gestation.day)
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
save(rf.model, file = "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.1.output/2.1.rf.model.20251010.Robj")


var_imp <- varImp(rf.model, scale = TRUE)
print(var_imp)
plot(var_imp)

imp_df <- var_imp$importance
imp_df$Variable <- rownames(imp_df)

imp_df_sorted <- imp_df[order(imp_df$Overall, decreasing = T), ]
imp_df_sorted.plot <- imp_df_sorted[1:10,]
imp_df_sorted.plot <- imp_df_sorted.plot[order(imp_df_sorted.plot$Overall, decreasing = F), ]
imp_df_sorted.plot$Variable <- factor(imp_df_sorted.plot$Variable, levels = imp_df_sorted.plot$Variable)
p <- ggplot(data=imp_df_sorted.plot)+
  geom_bar(stat = 'identity', mapping=aes(y = Variable, x = Overall), fill = "#af93c4")+
  scale_x_continuous(expand = c(0,0), limits=c(0,100))+
  theme(axis.text.x = element_text(colour = "black", angle = 30, hjust = 1), 
        axis.text.y = element_text(colour = "black", hjust = 1),
        legend.position = "none", panel.background = element_rect(fill = 'white'), 
        axis.line.x = element_line(colour = "black",  size=0.5, lineend = "butt"), 
        axis.line.y = element_line(colour = "black", size=0.5)) 
ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.1.output/2.1.importance.10.variable.rf.select.20251007.pdf",width = 8,height = 8)
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
save(pls_model, file = "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.1.output/2.1.pls_model.rf.select.20251010.Robj")


#####plot

best_ncomp <- pls_model$bestTune$ncomp
pred_train <- pls_model$pred
pred_train_best <- pred_train %>%
  dplyr::filter(ncomp == best_ncomp)

train_pred_df <- train.data %>%
  dplyr::mutate(
    Predicted_gestation.day = pred_train_best$pred[match(
      1:nrow(train.data),
      pred_train_best$rowIndex
    )]
  )
train_pred_df <- train_pred_df[,c("gestation.day" ,         
                                  "Predicted_gestation.day")]
colnames(train_pred_df) <- c("Chronological.age", "Predicted.age")
#######
p <- ggscatter(train_pred_df, x = "Chronological.age", y = "Predicted.age",          
               add = "reg.line",           
               conf.int = TRUE,    
               conf.int.level = 0.99,
               add.params = list(color = "#254689",                             
                                 fill = "#a0c9e5"),         
               cor.coef = TRUE,
               cor.coeff.args = list(method = "pearson", label.x = 120, label.sep = "\n"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()  
  )+ggtitle("Discovery Cohort (PLS model)")
ggsave (p,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.1.output/2.1.ggscatter.pls.model.rf.select.variant10.train.data.by.pls.20251010.pdf",width = 8,height = 8)

##################################################
#######test#######################################
##################################################
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/2.0.output/2.0.amniocytes.test.with.celltype.cor.0.8.exp.0.1.df.20251007.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/0.0.output/amnio.test.Robj")

gestation.day.test <- amniocytes.test@meta.data
gestation.day.test <- gestation.day.test[,c("gestation.day", "orig.ident")]
gestation.day.test <- gestation.day.test %>%
  group_by(orig.ident) %>%
  summarise(
    gestation.day = dplyr::first(gestation.day)
  )

df.test <- merge(df.test,gestation.day.test,by="orig.ident")
##################################################
test.data <- df.test
test.data[is.na(test.data)] = 0
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

test.data <- test.data[,c("gestation.day", "Predicted_Gestation.day")]
colnames(test.data) <- c("Chronological.age", "Predicted.age")
#######
p <- ggscatter(test.data, x = "Chronological.age", y = "Predicted.age",          
               add = "reg.line",           
               conf.int = TRUE,    
               conf.int.level = 0.99,
               add.params = list(color = "#B4355E",                             
                                 fill = "#f3dee0"),       
               cor.coef = TRUE,
               cor.coeff.args = list(method = "pearson", label.x = 140, label.sep = "\n"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()  
  )+ggtitle("Validation Cohort (PLS model)")
ggsave (p,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.1.output/2.1.ggscatter.pls.model.rf.select.variant10.test.data.by.pls.20251010.pdf",width = 8,height = 8)
##################################################
#######df.fgr.cohort1#############################
##################################################
load("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/5.0.output/amniocytes.cohort1.FGR.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/2.0.output/2.0.amniocytes.FGR.cohort1.with.celltype.cor.0.8.exp.0.1.df.20251007.Robj")
##################
gestation.day.FGR.cor1 <- amniocytes.FGR.cohort1@meta.data
gestation.day.FGR.cor1 <- gestation.day.FGR.cor1[,c("gestation.day", "orig.ident")]
gestation.day.FGR.cor1 <- gestation.day.FGR.cor1 %>%
  group_by(orig.ident) %>%
  summarise(
    gestation.day = dplyr::first(gestation.day)
  )

df.FGR.cohort1 <- merge(df.FGR.cohort1,gestation.day.FGR.cor1,by="orig.ident")
###############
FGR.cohort1.data <- df.FGR.cohort1
FGR.cohort1.data[is.na(FGR.cohort1.data)] = 0
rownames(FGR.cohort1.data) <- FGR.cohort1.data$orig.ident
FGR.cohort1.data <- FGR.cohort1.data[,-1]
colnames(FGR.cohort1.data) <- gsub("-",".",colnames(FGR.cohort1.data))

pred_test <- predict(
  pls_model,
  newdata = FGR.cohort1.data
)
FGR.cohort1.data <- FGR.cohort1.data %>%
  dplyr::mutate(
    Predicted_Gestation.day = as.numeric(pred_test)
  )


##########################################################plot
train_pred_df$group <- "Normal"

FGR.cohort1.data <- FGR.cohort1.data[,c("gestation.day","Predicted_Gestation.day")]
colnames(FGR.cohort1.data) <- c("Chronological.age", "Predicted.age")
FGR.cohort1.data$group <- "FGR"

matrix <- rbind(train_pred_df,FGR.cohort1.data)
matrix$delta <- matrix$Predicted.age-matrix$Chronological.age

p <- ggscatter(data = subset(matrix, group == "Normal"), x = "Chronological.age", y = "Predicted.age",          
               add = "reg.line",           
               conf.int = TRUE,    
               conf.int.level = 0.99,
               add.params = list(color = "#254689",                             
                                 fill = "#a0c9e5"),      
               cor.coeff.args = list(method = "pearson", label.x = 120, label.sep = "\n"))+
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank()   
  )+ggtitle("Train data with FGR cohort1 (Pls model)")+
  geom_point(data = subset(matrix, group == "FGR"), aes(x = Chronological.age, y = Predicted.age), 
             color = "#8b88b7", size = 4, shape = 17)+
  geom_smooth(data = subset(matrix, group == "FGR"), 
              aes(x = Chronological.age, y = Predicted.age),
              method = "lm", 
              se = TRUE, 
              level = 0.99,
              color = "#8b88b7", 
              fill = "#c0bfd8") +
  geom_vline(xintercept = 195, color = "black", linetype = "dashed") 
ggsave (p,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.1.output/2.1.ggscatter.pls.model.rf.select.cor.0.8.exp.0.1.fgr.cohort1.data.with.variant.10.by.pls.20251010.pdf",width = 8,height = 8)

####################
p <- ggplot(matrix, aes(x=delta))+
  geom_density(data = subset(matrix, group == "Normal"), fill = "#a0c9e5",alpha = 0.8,size = 0) +    
  geom_density(data = subset(matrix, group == "FGR"), fill = "#8b88b7", alpha = 0.6,size = 0) +  
  labs(
    title = "Normal cohort1 VS FGR cohort1",
    x = "Delta",
    y = "Density"
  ) +theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),    
    axis.line = element_blank(),     
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 11)
  )+geom_vline(
    xintercept = -14,
    linetype = "dashed",   
    color = "black",
    linewidth = 0.8
  )
ggsave (p,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.1.output/2.1.density.pls.model.cor.0.8.exp.0.1.fgr.cohort1.data.by.pls.20251010.pdf",width = 8,height = 8)

###
matrix$label <- ifelse(matrix$group == "Normal", 0, 1)

roc_delta <- roc(
  response = matrix$label,
  predictor = matrix$delta,
  quiet = TRUE
)

plot(
  roc_delta,
  col = "steelblue",
  lwd = 2,
  main = "ROC curve Test Set (Normal cohort1 and FGR cohort1)"
)


##################################################
#######df.fgr.cohort2#############################
##################################################
load("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/5.0.output/amniocytes.cohort2.FGR.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/2.0.output/2.0.amniocytes.FGR.cohort2.with.celltype.cor.0.8.exp.0.1.df.20251007.Robj")
##################
gestation.day.FGR.cor2 <- amniocytes.FGR.cohort2@meta.data
gestation.day.FGR.cor2 <- gestation.day.FGR.cor2[,c("gestation.day", "orig.ident")]
gestation.day.FGR.cor2 <- gestation.day.FGR.cor2 %>%
  group_by(orig.ident) %>%
  summarise(
    gestation.day = dplyr::first(gestation.day)
  )

df.FGR.cohort2 <- merge(df.FGR.cohort2,gestation.day.FGR.cor2,by="orig.ident")
###############
FGR.cohort2.data <- df.FGR.cohort2
FGR.cohort2.data[is.na(FGR.cohort2.data)] = 0
rownames(FGR.cohort2.data) <- FGR.cohort2.data$orig.ident
FGR.cohort2.data <- FGR.cohort2.data[,-1]
colnames(FGR.cohort2.data) <- gsub("-",".",colnames(FGR.cohort2.data))

pred_fgr.cor2 <- predict(
  pls_model,
  newdata = FGR.cohort2.data
)
FGR.cohort2.data <- FGR.cohort2.data %>%
  dplyr::mutate(
    Predicted_Gestation.day = as.numeric(pred_fgr.cor2)
  )
cor(FGR.cohort2.data$gestation.day, FGR.cohort2.data$Predicted_Gestation.day)

FGR.cohort2.data <- FGR.cohort2.data[,c("gestation.day", "Predicted_Gestation.day")]
colnames(FGR.cohort2.data) <- c("Chronological.age", "Predicted.age")
FGR.cohort2.data$group <- "FGR"
test.data$group <- "Normal"
matrix.2 <- rbind(test.data, FGR.cohort2.data)
matrix.2$delta <- matrix.2$Predicted.age-matrix.2$Chronological.age

p <- ggscatter(data = subset(matrix.2, group == "Normal"), x = "Chronological.age", y = "Predicted.age",          
               add = "reg.line",           
               conf.int = TRUE,    
               conf.int.level = 0.99,
               add.params = list(color = "#B4355E",                             
                                 fill = "#f3dee0"),    
               cor.coeff.args = list(method = "pearson", label.x = 120, label.sep = "\n"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()   
  )+ggtitle("Test data with FGR cohort2 (Pls model)")+
  geom_point(data = subset(matrix.2, group == "FGR"), aes(x = Chronological.age, y = Predicted.age), 
             color = "#8b88b7", size = 4, shape = 17)+
  geom_smooth(data = subset(matrix.2, group == "FGR"), 
              aes(x = Chronological.age, y = Predicted.age),
              method = "lm", 
              se = TRUE, 
              level = 0.99,
              color = "#8b88b7", 
              fill = "#c0bfd8") +
  geom_vline(xintercept = 195, color = "black", linetype = "dashed") 
ggsave (p,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.1.output/2.1.ggscatter.pls.model.rf.select.cor.0.8.exp.0.1.fgr.cohort2.data.with.variant.10.by.pls.20251010.pdf",width = 8,height = 8)

p <- ggplot(matrix.2, aes(x=delta))+
  geom_density(data = subset(matrix.2, group == "Normal"), fill = "#f3dee0",alpha = 0.8,size = 0) +    
  geom_density(data = subset(matrix.2, group == "FGR"), fill = "#8b88b7", alpha = 0.6,size = 0) +  
  labs(
    title = "Normal cohort2 VS FGR cohort2",
    x = "Delta",
    y = "Density"
  ) +theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),    
    axis.line = element_blank(),     
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 11)
  )+geom_vline(
    xintercept = -14,
    linetype = "dashed",  
    color = "black",
    linewidth = 0.8
  )
ggsave (p,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.1.output/2.1.density.pls.model.cor.0.8.exp.0.1.fgr.cohort2.data.by.pls.20251010.pdf",width = 8,height = 8)
####
matrix.2$label <- ifelse(matrix.2$group == "Normal", 0, 1)

roc_delta <- roc(
  response = matrix.2$label,
  predictor = matrix.2$delta,  # 注意：delta 越小越像 FGR，所以取负
  quiet = TRUE
)

auc(roc_delta)
#0.8865
ci.auc(roc_delta)
#95% CI: 0.7621-1 (DeLong)
plot(
  roc_delta,
  col = "brown",
  lwd = 2,
  main = "ROC curve Test Set (Normal cohort2 and FGR cohort2)"
)

