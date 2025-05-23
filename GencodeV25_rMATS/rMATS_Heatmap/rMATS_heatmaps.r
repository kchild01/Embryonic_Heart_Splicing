library(ComplexHeatmap)
library(circlize)
library(dplyr)
set.seed(1993)

#set path
directory <- "~/"

setwd(directory)

#SE
mydf <- read.table('~/cotney_lab/heart/rMATS_heatmap/final_SE.MATS.JC_PSI.txt', header = T, sep = '\t', row.names = 1, quote = "", stringsAsFactor = F)

first <- mydf[, 1:3]
mean_first <- rowMeans(first, na.rm = TRUE)
mean_first_matrix <- matrix(mean_first, nrow = nrow(mydf), ncol = 1)
rownames(mean_first_matrix) <- rownames(mydf)

second <- mydf[, 4:6]
mean_second <- rowMeans(second, na.rm = TRUE)
mean_second_matrix <- matrix(mean_second, nrow = nrow(mydf), ncol = 1)
rownames(mean_second_matrix) <- rownames(mydf)

third <- mydf[, 7:9]
mean_third <- rowMeans(third, na.rm = TRUE)
mean_third_matrix <- matrix(mean_third, nrow = nrow(mydf), ncol = 1)
rownames(mean_third_matrix) <- rownames(mydf)

fourth <- mydf[, 10:12]
mean_fourth <- rowMeans(fourth, na.rm = TRUE)
mean_fourth_matrix <- matrix(mean_fourth, nrow = nrow(mydf), ncol = 1)
rownames(mean_fourth_matrix) <- rownames(mydf)

fifth <- mydf[, 13:15]
mean_fifth <- rowMeans(fifth, na.rm = TRUE)
mean_fifth_matrix <- matrix(mean_fifth, nrow = nrow(mydf), ncol = 1)
rownames(mean_fifth_matrix) <- rownames(mydf)

sixth <- mydf[, 16:18]
mean_sixth <- rowMeans(sixth, na.rm = TRUE)
mean_sixth_matrix <- matrix(mean_sixth, nrow = nrow(mydf), ncol = 1)
rownames(mean_sixth_matrix) <- rownames(mydf)

seventh <- mydf[, 19:21]
mean_seventh <- rowMeans(seventh, na.rm = TRUE)
mean_seventh_matrix <- matrix(mean_seventh, nrow = nrow(mydf), ncol = 1)
rownames(mean_seventh_matrix) <- rownames(mydf)

eighth <- mydf[, 22:24]
mean_eighth <- rowMeans(eighth, na.rm = TRUE)
mean_eighth_matrix <- matrix(mean_eighth, nrow = nrow(mydf), ncol = 1)
rownames(mean_eighth_matrix) <- rownames(mydf)

df <- cbind(mean_first_matrix, mean_second_matrix, mean_third_matrix, mean_fourth_matrix, mean_fifth_matrix, mean_sixth_matrix, mean_seventh_matrix, mean_eighth_matrix)

mat <- t(apply(df, 1, scale))
colnames(mat) = c("CS13", "CS16", "CS17", "CS18", "CS19", "CS20", "CS21", "CS23")

ht_heart = Heatmap(mat,
                   show_row_names = FALSE,
                   cluster_columns = F,
                   km = 3,
                   row_km_repeats = 100,
                   gap = unit(4, "mm"),
                   heatmap_legend_param = list(
                     title = "Heart PSI Z-Score",
                     legend_height = unit(4, "cm"),
                     title_position = "leftcenter-rot"),
                   column_title = 'Embryonic Heart Development'
)

pdf("~/cotney_lab/heart/rMATS_heatmap/3cluster/plot/rMATS_heatmap.pdf", width = 10, height = 8)
ht = draw(ht_heart)
dev.off() 

#Cluster IDs
r.dend <- row_dend(ht)  #If needed, extract row dendrogram
rcl.list <- row_order(ht)  #Extract clusters (output is a list)

lapply(rcl.list, function(x) length(x))  #check/confirm size gene clusters

library(magrittr) # needed to load the pipe function '%%'

clu_df <- lapply(names(rcl.list), function(i){
  out <- data.frame(spliceID = rownames(mat[rcl.list[[i]],]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
  do.call(rbind, .)


write.table(clu_df, file= "~/cotney_lab/heart/rMATS_heatmap/3cluster/splice_clusters3.txt", sep="\t", quote=F, row.names=FALSE)

#A3SS
mydf <- read.table('./cotney_lab/heart/rMATS_heatmap/A3SS/final_A3SS.MATS.JC_PSI.txt', header = T, sep = '\t', row.names = 1, quote = "", stringsAsFactor = F)

first <- mydf[, 1:3]
mean_first <- rowMeans(first, na.rm = TRUE)
mean_first_matrix <- matrix(mean_first, nrow = nrow(mydf), ncol = 1)
rownames(mean_first_matrix) <- rownames(mydf)

second <- mydf[, 4:6]
mean_second <- rowMeans(second, na.rm = TRUE)
mean_second_matrix <- matrix(mean_second, nrow = nrow(mydf), ncol = 1)
rownames(mean_second_matrix) <- rownames(mydf)

third <- mydf[, 7:9]
mean_third <- rowMeans(third, na.rm = TRUE)
mean_third_matrix <- matrix(mean_third, nrow = nrow(mydf), ncol = 1)
rownames(mean_third_matrix) <- rownames(mydf)

fourth <- mydf[, 10:12]
mean_fourth <- rowMeans(fourth, na.rm = TRUE)
mean_fourth_matrix <- matrix(mean_fourth, nrow = nrow(mydf), ncol = 1)
rownames(mean_fourth_matrix) <- rownames(mydf)

fifth <- mydf[, 13:15]
mean_fifth <- rowMeans(fifth, na.rm = TRUE)
mean_fifth_matrix <- matrix(mean_fifth, nrow = nrow(mydf), ncol = 1)
rownames(mean_fifth_matrix) <- rownames(mydf)

sixth <- mydf[, 16:18]
mean_sixth <- rowMeans(sixth, na.rm = TRUE)
mean_sixth_matrix <- matrix(mean_sixth, nrow = nrow(mydf), ncol = 1)
rownames(mean_sixth_matrix) <- rownames(mydf)

seventh <- mydf[, 19:21]
mean_seventh <- rowMeans(seventh, na.rm = TRUE)
mean_seventh_matrix <- matrix(mean_seventh, nrow = nrow(mydf), ncol = 1)
rownames(mean_seventh_matrix) <- rownames(mydf)

eighth <- mydf[, 22:24]
mean_eighth <- rowMeans(eighth, na.rm = TRUE)
mean_eighth_matrix <- matrix(mean_eighth, nrow = nrow(mydf), ncol = 1)
rownames(mean_eighth_matrix) <- rownames(mydf)

df <- cbind(mean_first_matrix, mean_second_matrix, mean_third_matrix, mean_fourth_matrix, mean_fifth_matrix, mean_sixth_matrix, mean_seventh_matrix, mean_eighth_matrix)

mat <- t(apply(df, 1, scale))
colnames(mat) = c("CS13", "CS16", "CS17", "CS18", "CS19", "CS20", "CS21", "CS23")

ht_heart = Heatmap(mat,
                   show_row_names = FALSE,
                   cluster_columns = F,
                   km = 3,
                   row_km_repeats = 100,
                   gap = unit(4, "mm"),
                   heatmap_legend_param = list(
                     title = "Heart PSI Z-Score",
                     legend_height = unit(4, "cm"),
                     title_position = "leftcenter-rot"),
                   column_title = 'Embryonic Heart Development'
)

pdf("cotney_lab/heart/rMATS_heatmap/plot/A3SS/A3SS_rMATS_heatmap.pdf", width = 10, height = 8)
ht = draw(ht_heart)
dev.off() 

r.dend <- row_dend(ht)  #If needed, extract row dendrogram
rcl.list <- row_order(ht)  #Extract clusters (output is a list)

lapply(rcl.list, function(x) length(x))  #check/confirm size gene clusters

library(magrittr) # needed to load the pipe function '%%'

clu_df <- lapply(names(rcl.list), function(i){
  out <- data.frame(spliceID = rownames(mat[rcl.list[[i]],]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
  do.call(rbind, .)


write.table(clu_df, file= "cotney_lab/heart/rMATS_heatmap/A3SS/A3SS_splice_clusters.txt", sep="\t", quote=F, row.names=FALSE)

#A5SS
mydf <- read.table('./cotney_lab/heart/rMATS_heatmap/A5SS/final_A5SS.MATS.JC_PSI.txt', header = T, sep = '\t', row.names = 1, quote = "", stringsAsFactor = F)

first <- mydf[, 1:3]
mean_first <- rowMeans(first, na.rm = TRUE)
mean_first_matrix <- matrix(mean_first, nrow = nrow(mydf), ncol = 1)
rownames(mean_first_matrix) <- rownames(mydf)

second <- mydf[, 4:6]
mean_second <- rowMeans(second, na.rm = TRUE)
mean_second_matrix <- matrix(mean_second, nrow = nrow(mydf), ncol = 1)
rownames(mean_second_matrix) <- rownames(mydf)

third <- mydf[, 7:9]
mean_third <- rowMeans(third, na.rm = TRUE)
mean_third_matrix <- matrix(mean_third, nrow = nrow(mydf), ncol = 1)
rownames(mean_third_matrix) <- rownames(mydf)

fourth <- mydf[, 10:12]
mean_fourth <- rowMeans(fourth, na.rm = TRUE)
mean_fourth_matrix <- matrix(mean_fourth, nrow = nrow(mydf), ncol = 1)
rownames(mean_fourth_matrix) <- rownames(mydf)

fifth <- mydf[, 13:15]
mean_fifth <- rowMeans(fifth, na.rm = TRUE)
mean_fifth_matrix <- matrix(mean_fifth, nrow = nrow(mydf), ncol = 1)
rownames(mean_fifth_matrix) <- rownames(mydf)

sixth <- mydf[, 16:18]
mean_sixth <- rowMeans(sixth, na.rm = TRUE)
mean_sixth_matrix <- matrix(mean_sixth, nrow = nrow(mydf), ncol = 1)
rownames(mean_sixth_matrix) <- rownames(mydf)

seventh <- mydf[, 19:21]
mean_seventh <- rowMeans(seventh, na.rm = TRUE)
mean_seventh_matrix <- matrix(mean_seventh, nrow = nrow(mydf), ncol = 1)
rownames(mean_seventh_matrix) <- rownames(mydf)

eighth <- mydf[, 22:24]
mean_eighth <- rowMeans(eighth, na.rm = TRUE)
mean_eighth_matrix <- matrix(mean_eighth, nrow = nrow(mydf), ncol = 1)
rownames(mean_eighth_matrix) <- rownames(mydf)

df <- cbind(mean_first_matrix, mean_second_matrix, mean_third_matrix, mean_fourth_matrix, mean_fifth_matrix, mean_sixth_matrix, mean_seventh_matrix, mean_eighth_matrix)

mat <- t(apply(df, 1, scale))
colnames(mat) = c("CS13", "CS16", "CS17", "CS18", "CS19", "CS20", "CS21", "CS23")

ht_heart = Heatmap(mat,
                   show_row_names = FALSE,
                   cluster_columns = F,
                   km = 3,
                   row_km_repeats = 100,
                   gap = unit(4, "mm"),
                   heatmap_legend_param = list(
                     title = "Heart PSI Z-Score",
                     legend_height = unit(4, "cm"),
                     title_position = "leftcenter-rot"),
                   column_title = 'Embryonic Heart Development'
)

pdf("cotney_lab/heart/rMATS_heatmap/plot/A5SS/A5SS_rMATS_heatmap.pdf", width = 10, height = 8)
ht = draw(ht_heart)
dev.off() 

r.dend <- row_dend(ht)  #If needed, extract row dendrogram
rcl.list <- row_order(ht)  #Extract clusters (output is a list)

lapply(rcl.list, function(x) length(x))  #check/confirm size gene clusters

library(magrittr) # needed to load the pipe function '%%'

clu_df <- lapply(names(rcl.list), function(i){
  out <- data.frame(spliceID = rownames(mat[rcl.list[[i]],]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
  do.call(rbind, .)


write.table(clu_df, file= "cotney_lab/heart/rMATS_heatmap/A5SS/A5SS_splice_clusters.txt", sep="\t", quote=F, row.names=FALSE)

#MXE
mydf <- read.table('./cotney_lab/heart/rMATS_heatmap/MXE/final_MXE.MATS.JC_PSI.txt', header = T, sep = '\t', row.names = 1, quote = "", stringsAsFactor = F)

first <- mydf[, 1:3]
mean_first <- rowMeans(first, na.rm = TRUE)
mean_first_matrix <- matrix(mean_first, nrow = nrow(mydf), ncol = 1)
rownames(mean_first_matrix) <- rownames(mydf)

second <- mydf[, 4:6]
mean_second <- rowMeans(second, na.rm = TRUE)
mean_second_matrix <- matrix(mean_second, nrow = nrow(mydf), ncol = 1)
rownames(mean_second_matrix) <- rownames(mydf)

third <- mydf[, 7:9]
mean_third <- rowMeans(third, na.rm = TRUE)
mean_third_matrix <- matrix(mean_third, nrow = nrow(mydf), ncol = 1)
rownames(mean_third_matrix) <- rownames(mydf)

fourth <- mydf[, 10:12]
mean_fourth <- rowMeans(fourth, na.rm = TRUE)
mean_fourth_matrix <- matrix(mean_fourth, nrow = nrow(mydf), ncol = 1)
rownames(mean_fourth_matrix) <- rownames(mydf)

fifth <- mydf[, 13:15]
mean_fifth <- rowMeans(fifth, na.rm = TRUE)
mean_fifth_matrix <- matrix(mean_fifth, nrow = nrow(mydf), ncol = 1)
rownames(mean_fifth_matrix) <- rownames(mydf)

sixth <- mydf[, 16:18]
mean_sixth <- rowMeans(sixth, na.rm = TRUE)
mean_sixth_matrix <- matrix(mean_sixth, nrow = nrow(mydf), ncol = 1)
rownames(mean_sixth_matrix) <- rownames(mydf)

seventh <- mydf[, 19:21]
mean_seventh <- rowMeans(seventh, na.rm = TRUE)
mean_seventh_matrix <- matrix(mean_seventh, nrow = nrow(mydf), ncol = 1)
rownames(mean_seventh_matrix) <- rownames(mydf)

eighth <- mydf[, 22:24]
mean_eighth <- rowMeans(eighth, na.rm = TRUE)
mean_eighth_matrix <- matrix(mean_eighth, nrow = nrow(mydf), ncol = 1)
rownames(mean_eighth_matrix) <- rownames(mydf)

df <- cbind(mean_first_matrix, mean_second_matrix, mean_third_matrix, mean_fourth_matrix, mean_fifth_matrix, mean_sixth_matrix, mean_seventh_matrix, mean_eighth_matrix)

mat <- t(apply(df, 1, scale))
colnames(mat) = c("CS13", "CS16", "CS17", "CS18", "CS19", "CS20", "CS21", "CS23")

ht_heart = Heatmap(mat,
                   show_row_names = FALSE,
                   cluster_columns = F,
                   km = 3,
                   row_km_repeats = 100,
                   gap = unit(4, "mm"),
                   heatmap_legend_param = list(
                     title = "Heart PSI Z-Score",
                     legend_height = unit(4, "cm"),
                     title_position = "leftcenter-rot"),
                   column_title = 'Embryonic Heart Development'
)

pdf("cotney_lab/heart/rMATS_heatmap/plot/MXE/MXE_rMATS_heatmap.pdf", width = 10, height = 8)
ht = draw(ht_heart)
dev.off() 

r.dend <- row_dend(ht)  #If needed, extract row dendrogram
rcl.list <- row_order(ht)  #Extract clusters (output is a list)

lapply(rcl.list, function(x) length(x))  #check/confirm size gene clusters

library(magrittr) # needed to load the pipe function '%%'

clu_df <- lapply(names(rcl.list), function(i){
  out <- data.frame(spliceID = rownames(mat[rcl.list[[i]],]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
  do.call(rbind, .)


write.table(clu_df, file= "cotney_lab/heart/rMATS_heatmap/MXE/MXE_splice_clusters.txt", sep="\t", quote=F, row.names=FALSE)

#RI
mydf <- read.table('./cotney_lab/heart/rMATS_heatmap/RI/final_RI.MATS.JC_PSI.txt', header = T, sep = '\t', row.names = 1, quote = "", stringsAsFactor = F)

first <- mydf[, 1:3]
mean_first <- rowMeans(first, na.rm = TRUE)
mean_first_matrix <- matrix(mean_first, nrow = nrow(mydf), ncol = 1)
rownames(mean_first_matrix) <- rownames(mydf)

second <- mydf[, 4:6]
mean_second <- rowMeans(second, na.rm = TRUE)
mean_second_matrix <- matrix(mean_second, nrow = nrow(mydf), ncol = 1)
rownames(mean_second_matrix) <- rownames(mydf)

third <- mydf[, 7:9]
mean_third <- rowMeans(third, na.rm = TRUE)
mean_third_matrix <- matrix(mean_third, nrow = nrow(mydf), ncol = 1)
rownames(mean_third_matrix) <- rownames(mydf)

fourth <- mydf[, 10:12]
mean_fourth <- rowMeans(fourth, na.rm = TRUE)
mean_fourth_matrix <- matrix(mean_fourth, nrow = nrow(mydf), ncol = 1)
rownames(mean_fourth_matrix) <- rownames(mydf)

fifth <- mydf[, 13:15]
mean_fifth <- rowMeans(fifth, na.rm = TRUE)
mean_fifth_matrix <- matrix(mean_fifth, nrow = nrow(mydf), ncol = 1)
rownames(mean_fifth_matrix) <- rownames(mydf)

sixth <- mydf[, 16:18]
mean_sixth <- rowMeans(sixth, na.rm = TRUE)
mean_sixth_matrix <- matrix(mean_sixth, nrow = nrow(mydf), ncol = 1)
rownames(mean_sixth_matrix) <- rownames(mydf)

seventh <- mydf[, 19:21]
mean_seventh <- rowMeans(seventh, na.rm = TRUE)
mean_seventh_matrix <- matrix(mean_seventh, nrow = nrow(mydf), ncol = 1)
rownames(mean_seventh_matrix) <- rownames(mydf)

eighth <- mydf[, 22:24]
mean_eighth <- rowMeans(eighth, na.rm = TRUE)
mean_eighth_matrix <- matrix(mean_eighth, nrow = nrow(mydf), ncol = 1)
rownames(mean_eighth_matrix) <- rownames(mydf)

df <- cbind(mean_first_matrix, mean_second_matrix, mean_third_matrix, mean_fourth_matrix, mean_fifth_matrix, mean_sixth_matrix, mean_seventh_matrix, mean_eighth_matrix)

mat <- t(apply(df, 1, scale))
colnames(mat) = c("CS13", "CS16", "CS17", "CS18", "CS19", "CS20", "CS21", "CS23")

ht_heart = Heatmap(mat,
                   show_row_names = FALSE,
                   cluster_columns = F,
                   km = 3,
                   row_km_repeats = 100,
                   gap = unit(4, "mm"),
                   heatmap_legend_param = list(
                     title = "Heart PSI Z-Score",
                     legend_height = unit(4, "cm"),
                     title_position = "leftcenter-rot"),
                   column_title = 'Embryonic Heart Development'
)

pdf("cotney_lab/heart/rMATS_heatmap/plot/RI/RI_rMATS_heatmap.pdf", width = 10, height = 8)
ht = draw(ht_heart)
dev.off() 

r.dend <- row_dend(ht)  #If needed, extract row dendrogram
rcl.list <- row_order(ht)  #Extract clusters (output is a list)

lapply(rcl.list, function(x) length(x))  #check/confirm size gene clusters

library(magrittr) # needed to load the pipe function '%%'

clu_df <- lapply(names(rcl.list), function(i){
  out <- data.frame(spliceID = rownames(mat[rcl.list[[i]],]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
  do.call(rbind, .)


write.table(clu_df, file= "cotney_lab/heart/rMATS_heatmap/RI/RI_splice_clusters.txt", sep="\t", quote=F, row.names=FALSE)
