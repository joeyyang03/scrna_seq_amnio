rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/4.model/1.1.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 


library(pROC)
library(dplyr)
library(readr)
library(tidyverse)
library(ggrepel)
library(ggplot2)
library(reshape2)
library(caret)
library(randomForest)
library(stringr)
"%!in%" <- function(x,y)!('%in%'(x,y))
load("/share/home/yangjingyi/project/1.BD/up.load/4.model/1.0.output/1.0.amniocytes.normal.fgr.49.variant.df.20251007.Robj")
load("/share/home/yangjingyi/project/1.BD/up.load/4.model/1.0.output/1.0.amniocytes.49.variant.df.test.FGR.cohort2.20251009.Robj")
train.data <- df
rownames(train.data) <- train.data$orig.ident
train.data <- train.data %>% select(-orig.ident)
colnames(train.data) <- gsub("-","_",colnames(train.data))
test.data <- df.test
rownames(test.data) <- test.data$orig.ident
colnames(test.data) <- gsub("-","_",colnames(test.data))
train.data$group <- factor(train.data$group, levels = c("Normal", "FGR"))
levels(train.data$group)



#####################################################################################################

set.seed(123)
# ---- Step 1: RFE ----
ctrl <- rfeControl(functions = rfFuncs, method = "cv", number = 5)
rfe_model <- rfe(
  x = train.data[, -50],
  y = train.data$group,
  sizes = 1:49,
  rfeControl = ctrl
)
print(rfe_model)


rfe_res <- rfe_model$results
p <- ggplot(rfe_res, aes(x = Variables, y = Accuracy)) +
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
  labs(x = "Number of Features", y = "Accuracy")
ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/up.load/4.model/1.1.output/1.1.features.number.20251007.pdf",width = 8,height = 8)

best_features <- predictors(rfe_model)
print(best_features)

# ---- Step 2 ----
formula_best <- as.formula(
  paste("group ~", paste(best_features, collapse = " + "))
)

rf_best <- train(
  formula_best,
  data = train.data,
  method = "rf",
  trControl = trainControl(method = "cv", number = 5,
                           classProbs = TRUE, summaryFunction = twoClassSummary),
  metric = "ROC",
  importance = TRUE
)


var_imp <- varImp(rf_best, scale = TRUE)
print(var_imp)
plot(var_imp)

imp_df <- var_imp$importance


imp_df$Variable <- rownames(imp_df)


imp_df_sorted <- imp_df[order(imp_df$FGR), ]
imp_df_sorted.plot <- imp_df_sorted
imp_df_sorted.plot$Variable <- factor(imp_df_sorted.plot$Variable, levels = imp_df_sorted.plot$Variable)
p <- ggplot(data=imp_df_sorted.plot)+
  geom_bar(stat = 'identity', mapping=aes(y = Variable, x = FGR), fill = "#af93c4")+
  scale_x_continuous(expand = c(0,0), limits=c(0,100))+
  theme(axis.text.x = element_text(colour = "black", angle = 30, hjust = 1), 
        axis.text.y = element_text(colour = "black", hjust = 1),
        legend.position = "none", panel.background = element_rect(fill = 'white'), 
        axis.line.x = element_line(colour = "black",  size=0.5, lineend = "butt"), 
        axis.line.y = element_line(colour = "black", size=0.5)) 

ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/up.load/4.model/1.1.output/1.1.importance.48.variable.20251007.pdf",width = 8,height = 8)

imp_df_sorted.plot <- imp_df_sorted.plot[43:48,]
imp_df_sorted.plot$Variable <- factor(imp_df_sorted.plot$Variable, levels = imp_df_sorted.plot$Variable)
p <- ggplot(data=imp_df_sorted.plot)+
  geom_bar(stat = 'identity', mapping=aes(y = Variable, x = FGR), fill = "#af93c4")+
  coord_cartesian(xlim = c(70, 100)) +  
  theme(axis.text.x = element_text(colour = "black", angle = 30, hjust = 1), 
        axis.text.y = element_text(colour = "black", hjust = 1),
        legend.position = "none", panel.background = element_rect(fill = 'white'), 
        axis.line.x = element_line(colour = "black",  size=0.5, lineend = "butt"), 
        axis.line.y = element_line(colour = "black", size=0.5)) 

ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/up.load/4.model/1.1.output/1.1.importance.top6.variable.20251007.pdf",width = 8,height = 8)


top_vars <- rownames(var_imp$importance)[order(rowMeans(var_imp$importance), decreasing = TRUE)][1:6]
print(top_vars)
# ---- Step 3 ----
formula_top <- as.formula(
  paste("group ~", paste(top_vars, collapse = " + "))
)

# ---- Step 4 ----
set.seed(123)  
rf_top <- train(
  formula_top,
  data = train.data,
  method = "rf",
  trControl = trainControl(method = "cv", number = 5,
                           classProbs = TRUE,
                           summaryFunction = twoClassSummary,
                           savePredictions = TRUE),  
  metric = "ROC",
  importance = TRUE
)

print(rf_top)
# ---- Step 5: AUC ----

rf_top$results$ROC



roc_obj <- roc(response = rf_top$pred$obs,
               predictor = rf_top$pred$FGR,
               levels = c("Normal", "FGR")) 


plot(roc_obj, col = "steelblue", main = "ROC Curve on Train data")
auc(roc_obj)

ci.auc(roc_obj)


pred_test <- predict(rf_top, newdata = test.data, type = "prob")


prob_FGR <- pred_test$FGR

true_labels <- test.data$group

roc_test <- roc(response = true_labels,
                predictor = prob_FGR,
                levels = c("Normal", "FGR"))

plot(roc_test, col = "brown", main = "ROC Curve on Test Set")
auc(roc_test)

ci.auc(roc_test)


#####
pred_labels <- predict(rf_top, newdata = test.data, type = "raw")

result_df <- data.frame(
  Sample = rownames(test.data),  
  True_Label = true_labels,
  Pred_Label = pred_labels,
  Prob_FGR = prob_FGR
)

head(result_df)

###
information <- read.csv("/share/home/yangjingyi/project/1.BD/data/rawdata.20251001/df.summary.FGR.20251008.csv")
information <- information[information$vali.group == "cohort2",]
information <- information[,c("indicator", "gestation.week", "Sample.ID.FGR")]
colnames(information) <- c("indicator", "gestation.week", "Sample")

result_df <- merge(result_df,information, by = "Sample", all = T)
result_df$gestation.week[which(str_detect(result_df$Sample, "^GW"))]  <- sub("_.*", "", result_df$Sample[which(str_detect(result_df$Sample, "^GW"))])

result_df$indicator[which(str_detect(result_df$Sample, "^GW"))]  <- "Other.indicators"

result_df[result_df$indicator == "Other.indicators",]$indicator  <- "Normal"
df_summary <- result_df %>%
  group_by(gestation.week) %>%
  summarise(
    ultrasound_acc = mean(indicator == True_Label),
    model_acc = mean(Pred_Label == True_Label)
  )
df_long <- df_summary %>%
  pivot_longer(cols = c(ultrasound_acc, model_acc),
               names_to = "Source",
               values_to = "Accuracy")

ggplot(df_long, aes(x = gestation.week, y = Accuracy, color = Source, group = Source)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_classic() +
  labs(y = "Accuracy", color = "Source")

#####
df_long$week_num <- as.numeric(gsub("GW", "", df_long$gestation.week))
p <- ggplot(df_long, aes(x = week_num, y = Accuracy, color = Source)) +
  geom_smooth(se = FALSE, method = "loess") +
  theme_classic() +
  labs(x = "Gestation Week", y = "Accuracy", color = "Source") +
  scale_color_manual(
    values = c(
      "ultrasound_acc" = "#FAD797",
      "model_acc" = "#A888AF"      
    )
  )
ggsave (p ,file= "/share/home/yangjingyi/project/1.BD/up.load/4.model/1.1.output/1.1.acc.result.cohort2.loess.20251027.pdf",width = 8,height = 8)

