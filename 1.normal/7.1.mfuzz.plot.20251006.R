rm(list = ls())
setwd("/share/home/yangjingyi/project/1.BD/up.load/1.normal/7.1.output/")
.libPaths("/share/home/yangjingyi/Software/Rlib") 

library(dplyr)
library("Mfuzz")
library(Seurat)
library(tibble)
library(Cairo)
library(ggplot2)
library(pheatmap)
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
celltype <- basename(input_file)
celltype <- gsub(".mfuzz.20251005.RData", "", celltype)
celltype <- gsub("7.1.", "", celltype)

load(input_file)

#######################
row_max <- apply(group_mean, 1, max)
group_mean_filtered <- group_mean[row_max > 1, ]

##########################################

centroids <- mfuzz_cluster$centers

cluster_annotations <- data.frame(
  Cluster = rownames(centroids),
  Peak_Time = apply(centroids, 1, function(x) colnames(centroids)[which.max(x)]),
  Max_Expression = apply(centroids, 1, max)
)
print(cluster_annotations)

cluster_annotations <- cluster_annotations[order(cluster_annotations$Peak_Time, decreasing = F),]
cluster_annotations$Cluster <- factor(cluster_annotations$Cluster, levels = cluster_annotations$Cluster)
print(cluster_annotations)

membership_data <- as.data.frame(mfuzz_cluster$membership)
membership_data <- membership_data[rownames(membership_data) %in% rownames(group_mean_filtered),]

membership_data_annotations <- data.frame(
  Gene = rownames(membership_data),
  Cluster = apply(membership_data, 1, function(x) colnames(membership_data)[which.max(x)]),
  Max_Expression = apply(membership_data, 1, max)
)

unique(membership_data_annotations$Cluster)
membership_data_annotations$Cluster <- factor(membership_data_annotations$Cluster, levels = cluster_annotations$Cluster)
membership_data_annotations <- membership_data_annotations[order(membership_data_annotations$Cluster), ]
membership_data_annotations <- membership_data_annotations[membership_data_annotations$Max_Expression > 0.7,]

#################################################################
membership_data_annotations <- merge(
  membership_data_annotations,
  cluster_annotations[, c("Cluster", "Peak_Time")],
  by = "Cluster",
  all.x = TRUE
)
membership_data_annotations$cell.type <- celltype
table(membership_data_annotations$Peak_Time)

write.csv(x = membership_data_annotations, row.names = TRUE, file = paste0("7.1.mfuzz.gene.gestation.membership.annotations.",celltype,".20251006.csv"))


