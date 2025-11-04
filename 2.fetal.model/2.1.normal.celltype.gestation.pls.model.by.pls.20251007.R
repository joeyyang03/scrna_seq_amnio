rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/2.model/2.1.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 

library(dplyr)
library(readr)
library(Seurat)
library(pls)
library(ggplot2)
library(ggpubr)
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/0.0.output/amnio.train.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/2.0.output/2.0.amniocytes.with.celltype.cor.top30.exp.0.1.df.20251007.Robj")

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


train.ges.day <- train.data[,c("gestation.day"), drop = FALSE]
x <- model.matrix(gestation.day~.,train.data)[,-1]
y <- train.data$gestation.day
pls_model <- plsr(y ~ x, data = train.data, ncomp = 3,  validation = "CV", segments = 10) 

summary(pls_model)
predictions <- predict(pls_model, newdata = x)


predictions_vector <- predictions[, , 1] 


correlation_coefficient <- cor(predictions_vector, y)

print(correlation_coefficient)

coef(pls_model)

####
coef <- as.data.frame(coef(pls_model))
colnames(coef) <- c("Coefficient")
coef$factor <- rownames(coef)
coef <- coef[order(abs(coef$Coefficient), decreasing = F),]
coef$factor <- factor(coef$factor, levels = coef$factor)
top_coef <- coef[21:30, ]

p1 <- ggplot(data=top_coef)+
  geom_bar(stat = 'identity', mapping=aes(y = factor, x = abs(Coefficient)), fill = "#af93c4")+
  coord_cartesian(xlim = c(20, 60)) +  
  theme(axis.text.x = element_text(colour = "black", angle = 30, hjust = 1), 
        axis.text.y = element_text(colour = "black", hjust = 1),
        legend.position = "none", panel.background = element_rect(fill = 'white'), 
        axis.line.x = element_line(colour = "black",  size=0.5, lineend = "butt"), 
        axis.line.y = element_line(colour = "black", size=0.5)) 

ggsave (p1 ,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.1.output/2.1.pls.model.coef.top10.20251010.pdf",width = 8,height = 8)

####
pred_values <- predict(pls_model, x)[, 1, 3] 
results <- data.frame(
  orig.ident = rownames(train.data),
  Chronological.age = y,
  Predicted.age = pred_values,
  group = "Normal.train")

#######
p2 <- ggscatter(results, x = "Chronological.age", y = "Predicted.age",          
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
ggsave (p2,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.1.output/2.1.ggscatter.pls.model.cor.top30.exp.0.1.train.data.by.pls.20251010.pdf",width = 8,height = 8)
####################################################################################################
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/0.0.output/amnio.test.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/2.0.output/2.0.amniocytes.test.with.celltype.cor.top30.exp.0.1.df.20251007.Robj")
##################################################
gestation.day.test <- amniocytes.test@meta.data
gestation.day.test <- gestation.day.test[,c("gestation.day", "orig.ident")]
gestation.day.test <- gestation.day.test %>%
  group_by(orig.ident) %>%
  summarise(
    gestation.day = first(gestation.day)
  )

df.test <- merge(df.test,gestation.day.test,by="orig.ident")
##################################################
test.data <- df.test
test.data[is.na(test.data)] = 0
rownames(test.data) <- test.data$orig.ident
test.data <- test.data[,-1]

x.test <- model.matrix(gestation.day~.,test.data)[,-1]
test.predictions <- predict(pls_model, newdata = x.test)
test.predictions <- as.data.frame(test.predictions[,,1] )
colnames(test.predictions) <- c("Predicted.age")
test.predictions$orig.ident <- rownames(test.predictions)
test.predictions <- merge(test.predictions,gestation.day.test,by="orig.ident")
colnames(test.predictions) <- c("orig.ident", "Predicted.age", "Chronological.age")
#############################
#plot
p3 <- ggscatter(test.predictions, x = "Chronological.age", y = "Predicted.age",          
                add = "reg.line",           
                conf.int = TRUE,    
                conf.int.level = 0.99,
                add.params = list(color = "#B4355E",                             
                                  fill = "#f3dee0"),      
                cor.coef = TRUE,
                cor.coeff.args = list(method = "pearson", label.x = 120, label.sep = "\n"))+
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank()  
  )+ggtitle("Validation Cohort (PLS model)")
ggsave (p3,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.1.output/2.1.ggscatter.pls.model.cor.top30.exp.0.1.test.data.by.pls.20251010.pdf",width = 8,height = 8)

###############################
####df.fgr.cohort1#############
load("/share/home/yangjingyi/project/1.BD/up.load/3.FGR/5.0.output/amniocytes.cohort1.FGR.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/2.model/2.0.output/2.0.amniocytes.FGR.cohort1.with.celltype.cor.top30.exp.0.1.df.20251007.Robj")
##################################################
gestation.day.FGR.cor1 <- amniocytes.FGR.cohort1@meta.data
gestation.day.FGR.cor1 <- gestation.day.FGR.cor1[,c("gestation.day", "orig.ident")]
gestation.day.FGR.cor1 <- gestation.day.FGR.cor1 %>%
  group_by(orig.ident) %>%
  summarise(
    gestation.day = first(gestation.day)
  )

df.FGR.cohort1 <- merge(df.FGR.cohort1,gestation.day.FGR.cor1,by="orig.ident")
##################################################
FGR.cohort1.data <- df.FGR.cohort1
FGR.cohort1.data[is.na(FGR.cohort1.data)] = 0
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
matrix <- rbind(results,FGR.cohort1.predictions)
matrix$delta <- matrix$Chronological.age-matrix$Predicted.age

p4 <- ggscatter(data = subset(matrix, group == "Normal.train"), x = "Chronological.age", y = "Predicted.age",          
                add = "reg.line",           
                conf.int = TRUE,    
                conf.int.level = 0.99,
                add.params = list(color = "#254689",                             
                                  fill = "#a0c9e5"),
                cor.coeff.args = list(method = "pearson", label.x = 120, label.sep = "\n"))+
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank()   
  )+ggtitle("Train data with FGR cohort1 (Pls model)")+
  geom_point(data = subset(matrix, group == "FGR.cohort1"), aes(x = Chronological.age, y = Predicted.age), 
             color = "#B5382C", size = 4, shape = 17)+
  geom_smooth(data = subset(matrix, group == "FGR.cohort1"), 
              aes(x = Chronological.age, y = Predicted.age),
              method = "lm", 
              se = TRUE, 
              level = 0.99,
              color = "#B5382C", 
              fill = "#F4A7A7") +
  geom_vline(xintercept = 195, color = "black", linetype = "dashed") 
ggsave (p4,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.1.output/2.1.ggscatter.pls.model.cor.top30.exp.0.1.fgr.cohort1.data.with.cor.by.pls.20251010.pdf",width = 8,height = 8)

matrix.stageII <- matrix[matrix$Chronological.age < 196, ]
matrix.stageIII <- matrix[matrix$Chronological.age > 195, ]
p5 <- ggplot(matrix.stageII, aes(x=delta))+
  geom_density(data = subset(matrix.stageII, group == "Normal.train"), fill = "#254689",alpha = 0.5,size = 0) +    
  geom_density(data = subset(matrix.stageII, group == "FGR.cohort1"), fill = "#B5382C", alpha = 0.5,size = 0) +  
  labs(
    title = "Stage II",
    x = "Delta",
    y = "Density"
  ) +theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),    
    axis.line = element_blank(),     
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(2, "mm"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 11),
    plot.title = element_text(
      hjust = 0, face = "bold", size = 13, margin = margin(b = 8)
    ),
    plot.margin = margin(8, 8, 8, 8)
  )
ggsave (p5,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.1.output/2.1.density.pls.model.cor.top30.exp.0.1.fgr.cohort1.data.stageII.by.pls.20251010.pdf",width = 8,height = 8)


p5 <- ggplot(matrix.stageIII, aes(x=delta))+
  geom_density(data = subset(matrix.stageIII, group == "Normal.train"), fill = "#254689",alpha = 0.5,size = 0) +    
  geom_density(data = subset(matrix.stageIII, group == "FGR.cohort1"), fill = "#B5382C", alpha = 0.5,size = 0) +  
  labs(
    title = "Stage III",
    x = "Delta",
    y = "Density"
  ) +theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),    
    axis.line = element_blank(),     
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(2, "mm"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 11),
    plot.title = element_text(
      hjust = 0, face = "bold", size = 13, margin = margin(b = 8)
    ),
    plot.margin = margin(8, 8, 8, 8)
  )
ggsave (p5,file= "/share/home/yangjingyi/project/1.BD/up.load/2.model/2.1.output/2.1.density.pls.model.cor.top30.exp.0.1.fgr.cohort1.data.stageIII.by.pls.20251010.pdf",width = 8,height = 8)
