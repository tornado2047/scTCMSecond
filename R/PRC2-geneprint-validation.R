getwd()
Sys.time()


################## Part1-准备&输入 ##########################
rm(list = ls())
setwd("E:/zxf/workspace/gene_entropy/HSC")
obj <- readRDS("hsc.Rds")


rm(list = ls())
setwd("E:/zxf/workspace/gene_entropy/pbmc")
obj <- readRDS("pbmc3k.Rds")


rm(list = ls())
setwd("E:/zxf/workspace/gene_entropy/hPSCs")
obj <- readRDS("hPSCs_time.Rds")

rm(list = ls())
setwd("E:/zxf/workspace/gene_entropy/hPSCs/hPSCs result")
obj <- readRDS("hPSCs.Rds")

rm(list = ls())
setwd("E:/zxf/workspace/gene_entropy/Retinal Bipolar Neurons")
obj <- readRDS("RBN-cluster3.Rds")
obj <- readRDS("RBN-cluster13.Rds")
obj <- readRDS("RBN-cluster11.Rds")
obj <- readRDS("RBN-cluster2.Rds")
obj <- readRDS("RBN.Rds")

rm(list = ls())
setwd("E:/zxf/workspace/gene_entropy/hPSCnaive")
obj <- readRDS("hPSCs_naive.Rds")
obj <- readRDS("hPSCs_merge.Rds")


rm(list = ls())
setwd("E:/zxf/workspace/gene_entropy/hPSCcellfate")
obj <- readRDS("hPSCs_cellfate.Rds")
obj <- readRDS("hPSCs_cellfate_np.Rds")


rm(list = ls())
setwd("E:/zxf/workspace/gene_entropy/hPSC_PRC2_H9")
obj <- readRDS("hPSCs_PRC2.Rds")


rm(list = ls())
setwd("E:/zxf/workspace/gene_entropy/hPSC_H9_csv")
obj <- readRDS("hPSCs_csv.Rds")


rm(list = ls())
setwd("E:/zxf/workspace/gene_entropy/heart")
obj <- readRDS("heart.Rds")


rm(list = ls())
setwd("E:/zxf/workspace/gene_entropy/brain 10X")
obj <- readRDS("brain.Rds")


rm(list = ls())
setwd("E:/zxf/workspace/gene_entropy/breast cancer")
obj <- readRDS("breast.Rds")


rm(list = ls())
setwd("E:/zxf/workspace/gene_entropy/islets")
obj <- readRDS("islets.Rds")


rm(list = ls())
setwd("E:/zxf/workspace/gene_entropy/lung adenocarcinoma")
obj <- readRDS("lung.Rds")

######### Part2-需要加载的包  #######################################

.libPaths()
library(Rtsne)
library(tidyverse)
library(ggraph)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(patchwork)
library(igraph)
library(philentropy)
library(infotheo)
library(Rtsne)
library(corrplot)
library(scales)
library(philentropy)
library(reshape2)
library(stringr)
library(scatterplot3d)
library(factoextra)
library(plotly)
library(Seurat)

#########  Part3- 提取数据Extract statistically significant markers #######

obj.allmarkers <- FindAllMarkers(obj, min.pct = 0.2, only.pos = TRUE)
stat_marker <- subset(obj.allmarkers, abs(avg_log2FC) > 0.2)
# 去掉重复gene
stat_marker <- stat_marker[!duplicated(stat_marker$gene), ]

p <- DimPlot(obj, reduction = "umap")
h <- DimPlot(obj, reduction = "umap", 
             group.by = "orig.ident")
# cols = c("#225EA8","#984EA3"))

h|p


length(obj$seurat_clusters)
obj$seurat_clusters


### 函数1-2：获得自定义组合后的表达矩阵

get_custom_cluster_values <- function(obj, custom_clusters) {
  result_list <- list()
  for (new_cluster_num in 1:length(custom_clusters)) {
    # combine clusters based on custom_clusters
    combined_clusters <- custom_clusters[[new_cluster_num]]
    
    # filter stat_marker by combined clusters
    DEGs <- subset(stat_marker, cluster %in% combined_clusters)
    gene_names <- DEGs$gene
    marker <- data.frame(gene_names, DEGs$avg_log2FC)
    
    # get markers matrix for the combined clusters
    markers_matrix <- as.matrix(GetAssayData(obj, slot = "counts")[gene_names, WhichCells(obj, ident = combined_clusters)])
    
    M <- markers_matrix
    # log-scale 标准化数据
    A <- log2(M + 1)
    # 归一化数据
    A_norm <- t(apply(A, 1, function(x) x/sum(x)))
    
    R <- cor(t(A))
    
    # add the combined clusters' values to the result list
    result_list[[new_cluster_num]] <- list(DEGs = DEGs, marker = marker, M = M, A = A, A_norm =A_norm, R = R)
  }
  return(result_list)
}

custom_clusters <- list(2, c(0,1,3))

K_matstore <- get_custom_cluster_values(obj, custom_clusters)


length(K_matstore)


#### 函数2：计算sig——power
SigPowerFunction <- function(A_norm,R,DEGs) {
  # 1.Define the compute_kl_divergence function
  compute_kl_divergence <- function(matrix){
    n <- nrow(matrix)
    k <- ncol(matrix)
    uniform_mean <- apply(matrix, 1, mean) # Define uniform distribution as mean of each column
    
    S <- numeric(n) # Initialize S to a vector of zeros
    for (i in 1:n) {
      #归一化
      uniform_distribution <- as.vector(rep(uniform_mean[i], k))
      row_vec <- as.vector(matrix[i,])
      
      uniform_distribution <- uniform_distribution/sum(uniform_distribution)
      row_vec <-  row_vec/sum(row_vec)
      
      S[i] <- KL(rbind(row_vec, uniform_distribution))
    }
    S <- data.frame(S)
    row.names(S) <- row.names(matrix)
    return (S)
  }
  # 2.Compute KL divergence and store in S
  S <- compute_kl_divergence(A_norm)
  
  colnames(S) <- "S_value"
  
  # 3.计算每个基因的表达百分比,将结果存储为数据框
  percent_expressed <- rowSums(A_norm != 0) / ncol(A_norm) * 100
  Exp_percent  <- data.frame(Gene = rownames(A_norm), Exp_percent = percent_expressed)
  
  # 4. 整合数据框
  df <- data.frame(S,Exp_percent, DEGs)
  # Sort the data frame in descending order of S
  df_sorted <- df %>% arrange(desc(S))
  # Select the columns for the scatter plot
  df_part <- df_sorted[, c("avg_log2FC", "S_value","Exp_percent")]
  # Get the top 10 rows of the sorted data frame
  top_10_rows <- df_sorted %>% slice_head(n = 10)
  # 计算归一化后的 df_select$S_value
  # normalized_S_value <- (df_select$S_value - min(df_select$S_value)) / (max(df_select$S_value) - min(df_select$S_value))
  ratio_S_value <- df_part$S_value / sum(df_part$S_value)
  # 更新 df_select$MSI
  df_part$MSI <- df_part$avg_log2FC *ratio_S_value*1000
  # # Sort the data frame by the 'product' column in descending order
  # df_select_sig_power <- df_select %>% arrange(desc(sig_power))
  Sig_power_n <- df_part
  # 表达比例的threshold如何定？
  df_select <- df_part[df_part$Exp_percent > 50,]
  
  MSI <- subset(df_select,select="MSI")
  MSI <- MSI[order(MSI$MSI, decreasing = TRUE), , drop = FALSE]
  
  # top_genes <- head(rownames(MSI), 100)
  # 5.定义一个函数sparse_matrix_D，输入为向量S和常数k
  
  sparse_matrix_R <- function(S, k, R) {
    # 计算S中每个元素的k近邻
    dist_mat <- as.matrix(dist(S, diag = TRUE, upper = TRUE))
    knn_mat <- t(apply(dist_mat, 1, function(x) order(x)[1:k]))
    
    # 计算每个元素的LOF值
    lrd <- rep(0, length(S$MSI))
    length(lrd)
    for (i in 1:length(S$MSI)) {
      N_i <- knn_mat[i, ]
      dist_i <- dist_mat[i, N_i]
      lrd[i] <- k / sum(dist_i)
    }
    
    lof <- rep(0, length(S$MSI))
    for (i in 1:length(S$MSI)){
      N_i <- knn_mat[i, ]
      lrd_Ni <- lrd[N_i]
      dist_i <- dist_mat[i, N_i]
      lof[i] <- ifelse(sum(is.na(dist_i)) > 0, NA, sum(lrd_Ni / lrd[i] * dist_i) / k)
    }
    
    # 计算LOF值的阈值
    lof_threshold <- quantile(lof, 0.8, na.rm = TRUE)
    
    # 将S稀疏化
    S_sparse <- S$MSI * (lof >= lof_threshold)
    S_sparse <- data.frame(S_sparse)
    colnames(S_sparse) <-"MSI"
    rownames(S_sparse) <- rownames(S)
    S_sparse <- S_sparse[S_sparse != 0, , drop = FALSE]
    # 将 S_sparse 转换为逻辑向量
    important_elements <- rownames(R) %in% rownames(S_sparse)
    # 使用 R * important_elements 进行稀疏化
    R_sparse <- R * important_elements
    diag(R_sparse) <- 0
    R_sparse[R_sparse< 0.2] <- 0
    # result_list <- list(R=R_sparse,S_sparse = S_sparse, lof_threshold = lof_threshold,important_elements=important_elements)
    return(R_sparse)
    # return(result_list)
  }
  
  sparse_corr_mat <- sparse_matrix_R(MSI,10,R)
  
  create_msi_network <- function(R, y) {
    # Create the graph
    G <- graph.adjacency(R, mode = "undirected", weighted = TRUE)
    G <- delete.vertices(G, which(degree(G) == 0))
    
    # Set node sizes based on y values
    default_size <- 0.1
    for (node_name in V(G)$name) {
      row_index <- match(node_name, rownames(y))
      
      if (!is.na(row_index) && length(row_index) > 0) {
        V(G)$size[V(G)$name == node_name] <- as.numeric(y[row_index, ])
      } else {
        V(G)$size[V(G)$name == node_name] <- default_size
      }
    }
    
    # Set edge widths based on R values
    E(G)$width <- R[as.matrix(get.edgelist(G))]
    
    return(G)
  }
  
  # Call the function and check for debugging information
  G <- create_msi_network(sparse_corr_mat,MSI)
  
  # 7.将节点名称作为第一列输出
  deg <- degree(G)
  
  deg_df <- data.frame(degree = deg)
  
  df_select$Degree <- ifelse(row.names(df_select) %in% row.names(deg_df), deg_df[row.names(df_select), "degree"], 0)
  
  # colnames(df_select)[colnames(df_select)=="deg_df$degree"] <- "Degree"
  
  ####  得到重要数据框 @@@@@
  Sig_power <- df_select
  # Sig_power <- df_select[df_select$Exp_percent > 25,] ####  得到重要数据框 @@@@@
  
  # Sig_power <- df_select
  Sig_power_sorted <- Sig_power[order(-Sig_power$MSI),]
  
  top_names <- head(rownames(Sig_power_sorted), 20)
  
  # 6. create network of sparse_corr_mat via MSI
  result_list <- list(sparse_corr_mat = sparse_corr_mat,Sig_power = Sig_power, Sig_power_n = Sig_power_n,
                      MSI = MSI, ER_factor = top_names, G = G)
  return(result_list) 
  
}

### adjust function sigpower
SigPowerFunction <- function(A_norm,R,DEGs) {
  # 1.Define the compute_kl_divergence function
  compute_kl_divergence <- function(matrix){
    n <- nrow(matrix)
    k <- ncol(matrix)
    uniform_mean <- apply(matrix, 1, mean) # Define uniform distribution as mean of each column
    
    S <- numeric(n) # Initialize S to a vector of zeros
    for (i in 1:n) {
      #归一化
      uniform_distribution <- as.vector(rep(uniform_mean[i], k))
      row_vec <- as.vector(matrix[i,])
      
      uniform_distribution <- uniform_distribution/sum(uniform_distribution)
      row_vec <-  row_vec/sum(row_vec)
      
      S[i] <- KL(rbind(row_vec, uniform_distribution))
    }
    S <- data.frame(S)
    row.names(S) <- row.names(matrix)
    return (S)
  }
  # 2.Compute KL divergence and store in S
  S <- compute_kl_divergence(A_norm)
  
  colnames(S) <- "S_value"
  
  # 3.计算每个基因的表达百分比,将结果存储为数据框
  percent_expressed <- rowSums(A_norm != 0) / ncol(A_norm) * 100
  Exp_percent  <- data.frame(Gene = rownames(A_norm), Exp_percent = percent_expressed)
  
  # 4. 整合数据框
  df <- data.frame(S,Exp_percent, DEGs)
  # Sort the data frame in descending order of S
  df_sorted <- df %>% arrange(desc(S))
  # Select the columns for the scatter plot
  df_select <- df_sorted[, c("avg_log2FC", "S_value","Exp_percent")]
  # Get the top 10 rows of the sorted data frame
  top_10_rows <- df_sorted %>% slice_head(n = 10)
  # Add a new column to the data frame that calculates the product of 'avg_log2FC' and 'S'
  # df_select$sig_power <- df_select$avg_log2FC * ((log(df_select$degree,base = 2))^df_select$S)
  # 计算归一化后的 df_select$S_value
  # normalized_S_value <- (df_select$S_value - min(df_select$S_value)) / (max(df_select$S_value) - min(df_select$S_value))
  ratio_S_value <- df_select$S_value / sum(df_select$S_value)
  # 更新 df_select$MSI
  df_select$MSI <- df_select$avg_log2FC *ratio_S_value*1000
  # # Sort the data frame by the 'product' column in descending order
  # df_select_sig_power <- df_select %>% arrange(desc(sig_power))
  
  df_select <- df_select[df_select$Exp_percent > 50,]
  MSI <- subset(df_select,select="MSI")
  MSI <- MSI[order( MSI$MSI, decreasing = TRUE), , drop = FALSE]
  
  # top_genes <- head(rownames(MSI), 100)
  # 5.定义一个函数sparse_matrix_D，输入为向量S和常数k
  sparse_matrix_D <- function(S, k, R) {
    # 计算S中每个元素的k近邻
    dist_mat <- as.matrix(dist(S, diag = TRUE, upper = TRUE))
    knn_mat <- t(apply(dist_mat, 1, function(x) order(x)[1:k]))
    
    # 计算每个元素的LOF值
    lrd <- rep(0, length(S$MSI))
    length(lrd)
    for (i in 1:length(S$MSI)) {
      N_i <- knn_mat[i, ]
      dist_i <- dist_mat[i, N_i]
      lrd[i] <- k / sum(dist_i)
    }
    
    lof <- rep(0, length(S$MSI))
    for (i in 1:length(S$MSI)){
      N_i <- knn_mat[i, ]
      lrd_Ni <- lrd[N_i]
      dist_i <- dist_mat[i, N_i]
      lof[i] <- ifelse(sum(is.na(dist_i)) > 0, NA, sum(lrd_Ni / lrd[i] * dist_i) / k)
    }
    
    # 计算LOF值的阈值
    lof_threshold <- quantile(lof, 0.8, na.rm = TRUE)
    
    # 将S稀疏化
    S_sparse <- S$MSI * (lof >= lof_threshold)
    S_sparse <- data.frame(S_sparse)
    colnames(S_sparse) <-"MSI"
    
    # 自定义函数，用于计算差的绝对值
    abs_diff <- function(x, y){
      abs(x - y)
    }
    # 应用outer()函数来创建矩阵
    A <- outer(S_sparse$MSI, S_sparse$MSI, abs_diff)
    rownames(A) <- rownames(S_sparse)
    colnames(A) <- rownames(S_sparse)
    
    # 将每列元素小于0.95分位数的置为0
    A[A < 0.001] <- 0
    R[A==0] <- 0
    R[R < 0.3] <- 0
    diag(R) <- 0
    return(R)
    # return(list(S_sparse = S_sparse, lof_threshold = lof_threshold))
  }
  
  sparse_corr_mat <- sparse_matrix_D(MSI,10,R)
  
  create_gene_network <- function(R, y) {
    
    # Select the top 50 genes from y
    top_msi_genes <- rownames(head(y, 200))
    
    # Extract only the top 50 genes from y
    top_genes <- y[top_msi_genes, , drop = FALSE]
    
    # Print dimensions and structure for debugging
    print(dim(top_genes))
    print(str(top_genes))
    
    # Set non-top genes in R to 0
    R[!(rownames(R) %in% top_msi_genes), ] <- 0
    R[, !(colnames(R) %in% top_msi_genes)] <- 0
    
    # Print dimensions for debugging
    print(dim(R))
    
    # Create the graph
    G <- graph.adjacency(R, mode = "undirected", weighted = TRUE)
    G <- delete.vertices(G, which(degree(G) == 0))
    
    # Set vertex attributes
    for (i in 1:vcount(G)) {
      vertex_name <- V(G)$name[i]
      V(G)$size[i] <- as.numeric(top_genes$MSI[top_msi_genes == vertex_name])
    }
    
    # Print length of node_sizes for debugging
    print(length(V(G)$size))
    
    return(G)
  }
  
  # Call the function and check for debugging information
  G <- create_gene_network(sparse_corr_mat, MSI)
  # 7.将节点名称作为第一列输出
  deg <- degree(G)
  
  deg_df <- data.frame(degree = deg)
  
  df_select$Degree <- ifelse(row.names(df_select) %in% row.names(deg_df), deg_df[row.names(df_select), "degree"], 0)
  
  # colnames(df_select)[colnames(df_select)=="deg_df$degree"] <- "Degree"
  
  Sig_power <- df_select ####  得到重要数据框 @@@@@
  # Sig_power <- df_select[df_select$Exp_percent > 25,] ####  得到重要数据框 @@@@@
  
  # Sig_power <- df_select
  Sig_power_sorted <- Sig_power[order(-Sig_power$MSI),]
  
  top_names <- head(rownames(Sig_power_sorted), 20)
  
  # 6. create network of sparse_corr_mat via MSI
  result_list <- list(sparse_corr_mat = sparse_corr_mat,Sig_power = Sig_power, MSI = MSI, ER_factor = top_names, G = G)
  return(result_list) 
  
}

### 将S值设为归一化的信息量
SigPowerFunction <- function(A_norm,R,DEGs) {
  compute_normalized_information <- function(matrix) {
    n <- nrow(matrix)
    k <- ncol(matrix)
    
    # Calculate entropy for each row
    entropy <- numeric(n)
    for (i in 1:n) {
      row_vec <- as.vector(matrix[i, ])
      row_entropy <- -sum(row_vec * log(row_vec + 1e-10)) / log(k)
      entropy[i] <- row_entropy
    }
    
    # Normalize entropy values
    total_entropy <- sum(entropy)
    normalized_info <- entropy / total_entropy
    
    return(data.frame(S_value = normalized_info))
  }
  
  
  # 2.Compute KL divergence and store in S
  S <-compute_normalized_information(A_norm)
  
  colnames(S) <- "S_value"
  
  # 3.计算每个基因的表达百分比,将结果存储为数据框
  percent_expressed <- rowSums(A_norm != 0) / ncol(A_norm) * 100
  Exp_percent  <- data.frame(Gene = rownames(A_norm), Exp_percent = percent_expressed)
  
  # 4. 整合数据框
  df <- data.frame(S,Exp_percent, DEGs)
  # Sort the data frame in descending order of S
  df_sorted <- df %>% arrange(desc(S))
  # Select the columns for the scatter plot
  df_select <- df_sorted[, c("avg_log2FC", "S_value","Exp_percent")]
  # Get the top 10 rows of the sorted data frame
  top_10_rows <- df_sorted %>% slice_head(n = 10)
  # 更新 df_select$MSI
  df_select$MSI <- df_select$avg_log2FC *df_select$S_value*100
  # # Sort the data frame by the 'product' column in descending order
  # df_select_sig_power <- df_select %>% arrange(desc(sig_power))
  
  df_select <- df_select[df_select$Exp_percent > 20,]
  MSI <- subset(df_select,select="MSI")
  MSI <- MSI[order( MSI$MSI, decreasing = TRUE), , drop = FALSE]
  
  # top_genes <- head(rownames(MSI), 100)
  # 5.定义一个函数sparse_matrix_D，输入为向量S和常数k
  sparse_matrix_D <- function(S, k, R) {
    # 计算S中每个元素的k近邻
    dist_mat <- as.matrix(dist(S, diag = TRUE, upper = TRUE))
    knn_mat <- t(apply(dist_mat, 1, function(x) order(x)[1:k]))
    
    # 计算每个元素的LOF值
    lrd <- rep(0, length(S$MSI))
    length(lrd)
    for (i in 1:length(S$MSI)) {
      N_i <- knn_mat[i, ]
      dist_i <- dist_mat[i, N_i]
      lrd[i] <- k / sum(dist_i)
    }
    
    lof <- rep(0, length(S$MSI))
    for (i in 1:length(S$MSI)){
      N_i <- knn_mat[i, ]
      lrd_Ni <- lrd[N_i]
      dist_i <- dist_mat[i, N_i]
      lof[i] <- ifelse(sum(is.na(dist_i)) > 0, NA, sum(lrd_Ni / lrd[i] * dist_i) / k)
    }
    
    # 计算LOF值的阈值
    lof_threshold <- quantile(lof, 0.7, na.rm = TRUE)
    
    # 将S稀疏化
    S_sparse <- S$MSI * (lof >= lof_threshold)
    S_sparse <- data.frame(S_sparse)
    colnames(S_sparse) <-"MSI"
    
    # 自定义函数，用于计算差的绝对值
    abs_diff <- function(x, y){
      abs(x - y)
    }
    # 应用outer()函数来创建矩阵
    A <- outer(S_sparse$MSI, S_sparse$MSI, abs_diff)
    rownames(A) <- rownames(S_sparse)
    colnames(A) <- rownames(S_sparse)
    
    # 将每列元素小于0.95分位数的置为0
    A[A < 0.001] <- 0
    R[A==0] <- 0
    R[R < 0.5] <- 0
    diag(R) <- 0
    return(R)
    # return(list(S_sparse = S_sparse, lof_threshold = lof_threshold))
  }
  
  sparse_corr_mat <- sparse_matrix_D(MSI,10,R)
  
  create_gene_network <- function(R, y) {
    
    # Select the top 50 genes from y
    top_msi_genes <- rownames(head(y, 200))
    
    # Extract only the top 50 genes from y
    top_genes <- y[top_msi_genes, , drop = FALSE]
    
    # Print dimensions and structure for debugging
    print(dim(top_genes))
    print(str(top_genes))
    
    # Set non-top genes in R to 0
    R[!(rownames(R) %in% top_msi_genes), ] <- 0
    R[, !(colnames(R) %in% top_msi_genes)] <- 0
    
    # Print dimensions for debugging
    print(dim(R))
    
    # Create the graph
    G <- graph.adjacency(R, mode = "undirected", weighted = TRUE)
    G <- delete.vertices(G, which(degree(G) == 0))
    
    # Set vertex attributes
    for (i in 1:vcount(G)) {
      vertex_name <- V(G)$name[i]
      V(G)$size[i] <- as.numeric(top_genes$MSI[top_msi_genes == vertex_name])
    }
    
    # Print length of node_sizes for debugging
    print(length(V(G)$size))
    
    return(G)
  }
  
  # Call the function and check for debugging information
  G <- create_gene_network(sparse_corr_mat, MSI)
  # 7.将节点名称作为第一列输出
  deg <- degree(G)
  
  deg_df <- data.frame(degree = deg)
  
  df_select$Degree <- ifelse(row.names(df_select) %in% row.names(deg_df), deg_df[row.names(df_select), "degree"], 0)
  
  # colnames(df_select)[colnames(df_select)=="deg_df$degree"] <- "Degree"
  
  Sig_power <- df_select ####  得到重要数据框 @@@@@
  # Sig_power <- df_select[df_select$Exp_percent > 25,] ####  得到重要数据框 @@@@@
  
  # Sig_power <- df_select
  Sig_power_sorted <- Sig_power[order(-Sig_power$MSI),]
  
  top_names <- head(rownames(Sig_power_sorted), 30)
  
  # 6. create network of sparse_corr_mat via MSI
  result_list <- list(sparse_corr_mat = sparse_corr_mat,Sig_power = Sig_power, MSI = MSI, ER_factor = top_names, G = G)
  return(result_list) 
  
}




#################

library(infotheo)
library(entropy)
library(progress)

X <- K_matstore[[1]]$A

# 5.定义一个函数输入

### 计算表达矩阵的互信息
create_MI_mat <- function(expr_matrix) {
  pb <- progress_bar$new(format = "[:bar] :percent ETA: :eta", total = nrow(expr_matrix) * (nrow(expr_matrix) - 1) / 2)

  start_time <- Sys.time()  # 记录程序开始时间

  expr_matrix<- round(expr_matrix)
  n_genes <- nrow(expr_matrix)
  MI_matrix <- matrix(0, nrow = n_genes, ncol = n_genes)

  for (i in 1:(n_genes - 1)) {
    for (j in (i + 1):n_genes) {
      gene1_expr <- expr_matrix[i, ]
      gene2_expr <- expr_matrix[j, ]
      # 计算互信息
      MI <- mutinformation(gene1_expr,gene2_expr,method= "emp")

      MI_matrix[i, j] <- MI
      MI_matrix[j, i] <- MI
      pb$tick()  # 更新进度条
    }
  }

  diag(MI_matrix) <- 0

  end_time <- Sys.time()  # 记录程序结束时间
  elapsed_time <- end_time - start_time  # 计算程序运行时间

  cat("\n程序运行时间：", format(elapsed_time, units = "secs"), "\n")

  return(MI_matrix)
}


### 将S值设为归一化的信息量，使用互信息增益矩阵构造网络
SigPowerFunction2 <- function(MI_mat,DEGs) {

  compute_normalized_information <- function(matrix) {
    n <- nrow(matrix)
    k <- ncol(matrix)

    # Calculate entropy for each row
    entropy <- numeric(n)
    for (i in 1:n) {
      row_vec <- as.vector(matrix[i, ])
      row_entropy <- -sum(row_vec * log(row_vec + 1e-10)) / log(k)
      entropy[i] <- row_entropy
    }

    # Normalize entropy values
    total_entropy <- sum(entropy)
    normalized_info <- entropy / total_entropy

    return(data.frame(S_value = normalized_info))
  }

  # 2.Compute KL divergence and store in S
  S <- compute_normalized_information(A)
  colnames(S) <- "S_value"

  # 3.计算每个基因的表达百分比,将结果存储为数据框
  percent_expressed <- rowSums(A != 0) / ncol(A) * 100
  Exp_percent  <- data.frame(Gene = rownames(A), Exp_percent = percent_expressed)

  # 4. 整合数据框
  df <- data.frame(S,Exp_percent, DEGs)
  # Sort the data frame in descending order of S
  df_sorted <- df %>% arrange(desc(S))
  # Select the columns for the scatter plot
  df_select <- df_sorted[, c("avg_log2FC", "S_value","Exp_percent")]
  # Get the top 10 rows of the sorted data frame
  top_10_rows <- df_sorted %>% slice_head(n = 10)
  # 更新 df_select$MSI
  df_select$MSI <- df_select$avg_log2FC *df_select$S_value*100
  # # Sort the data frame by the 'product' column in descending order
  # df_select_sig_power <- df_select %>% arrange(desc(sig_power))
  df_select <- df_select[df_select$Exp_percent > 20,]
  MSI <- subset(df_select,select="MSI")
  MSI <- MSI[order( MSI$MSI, decreasing = TRUE), , drop = FALSE]

  create_gene_network <- function(MI_mat, y) {

    # Select the top 50 genes from y
    top_msi_genes <- rownames(head(y, 200))

    # Extract only the top 50 genes from y
    top_genes <- y[top_msi_genes, , drop = FALSE]

    # Print dimensions and structure for debugging
    print(dim(top_genes))
    print(str(top_genes))

    # Set non-top genes in R to 0
    MI_mat[!(rownames(MI_mat) %in% top_msi_genes), ] <- 0
    MI_mat[, !(colnames(MI_mat) %in% top_msi_genes)] <- 0

    # Print dimensions for debugging
    print(dim(MI_mat))

    # Create the graph
    G <- graph.adjacency(MI_mat, mode = "undirected", weighted = TRUE)
    G <- delete.vertices(G, which(degree(G) == 0))

    # Set vertex attributes
    for (i in 1:vcount(G)) {
      vertex_name <- V(G)$name[i]
      V(G)$size[i] <- as.numeric(top_genes$MSI[top_msi_genes == vertex_name])
    }

    # Print length of node_sizes for debugging
    print(length(V(G)$size))

    return(G)
  }

  # Call the function and check for debugging information
  G <- create_gene_network(MI_mat, MSI)
  # 7.将节点名称作为第一列输出
  deg <- degree(G)

  deg_df <- data.frame(degree = deg)

  df_select$Degree <- ifelse(row.names(df_select) %in% row.names(deg_df), deg_df[row.names(df_select), "degree"], 0)

  # colnames(df_select)[colnames(df_select)=="deg_df$degree"] <- "Degree"

  Sig_power <- df_select ####  得到重要数据框 @@@@@
  # Sig_power <- df_select[df_select$Exp_percent > 25,] ####  得到重要数据框 @@@@@

  # Sig_power <- df_select
  Sig_power_sorted <- Sig_power[order(-Sig_power$MSI),]

  top_names <- head(rownames(Sig_power_sorted), 30)

  # 6. create network of sparse_corr_mat via MSI
  result_list <- list(Sig_power = Sig_power, MSI = MSI, ER_factor = top_names, G = G)
  return(result_list)

}


getMI_results <- list()
for (i in 1:2) {
  MI_mat <- create_MI_mat(K_matstore[[i]]$A)
  getMI_results <- MI_mat
}


#################

# 绘制 sigpower plot
create_sig_power_plot <- function(Sig_power) {
  # sort the data by sig_power in descending order
  Sig_power_sorted1 <- Sig_power[order(-Sig_power$S_value),]
  Sig_power_sorted2 <- Sig_power[order(-Sig_power$avg_log2FC),]
  Sig_power_sorted3 <- Sig_power[order(-Sig_power$MSI),]

  # # Extract the top 20 names and corresponding sig_power values
  # top_names1 <- head(rownames(Sig_power_sorted1), 5)
  # top_names2 <- head(rownames(Sig_power_sorted2), 5)
  top_names3 <- head(rownames(Sig_power_sorted3), 20)


  # top_names <- c(top_names1, top_names2,  top_names3)
  # top_names <- unique(top_names)
  top_names <- top_names3

  # Merge the sorted data
  Sig_power_sorted <- rbind(Sig_power_sorted1, Sig_power_sorted2)

  # 设置字体主题
  my_theme <- theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold", family = "Times New Roman"),
      axis.text = element_text(face = "bold", family = "Times New Roman"),
      legend.title = element_text(face = "bold", family = "Times New Roman"),
      legend.text = element_text(face = "bold", family = "Times New Roman"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black")
    )
  Sig2DPlot <- ggplot(Sig_power_sorted, aes(x = avg_log2FC, y = MSI)) +
    geom_point(aes(size = Degree, fill = Exp_percent), shape = 21, color = "black") +
    scale_fill_gradient(low = "lightcyan3", high = "violetred") +
    scale_size(range = c(2, 10)) +
    labs(x = "avg_log2FC", y = "MSI", size = "Degree", fill = "Exp_percent" ) +
    my_theme +
    geom_text_repel(data = Sig_power_sorted[top_names, ],
                    aes(label = rownames(Sig_power_sorted[top_names, ])),
                    family = "Times New Roman", color = "black",
                    box.padding = 0.3, point.padding = 0.3, size = 3.5,
                    force = 0.5, segment.color = "black") +
    theme(legend.position = "top") +
    guides(fill = guide_legend(title = "Exp_percent",
                               key_height = unit(3, "cm"),
                               override.aes = list(size =5)))

  return(Sig2DPlot)
}









getsigpower_results <- list() 
for (i in 1:2) { 
  getsigpower <- SigPowerFunction2(K_matstore[[i]]$A, K_matstore[[i]]$DEGs) 
  getsigpower_results[[i]] <- getsigpower 
}


# 输出量Sig_power
Sig_power_list <- list() 
for (i in 1:2) { 
  Sig_power_list[[i]] <- getsigpower_results[[i]]$Sig_power
} 


# 输出量sig_power plot
sig_power_plot_list <- list() 
for (i in 1:2) { 
  sig_power_plot_list[[i]] <- create_sig_power_plot(Sig_power_list[[i]])
} 


df1 <- Sig_power_list[[1]]
df2 <- Sig_power_list[[2]]

sig_power_plot_list[[1]]
sig_power_plot_list[[2]]



######  分别输出三列的值


output_fc_top100 <- data.frame(avg_log2FC = numeric(0), cluster = numeric(0))
# 循环遍历列表元素
for (i in 1:length(getsigpower_results)) {
  Sig_power <- getsigpower_results[[i]]$Sig_power
  
  Sig_power <- Sig_power[order(-Sig_power$avg_log2FC),]
  fc_factor <- head(rownames(Sig_power), 100)
  cluster <- i - 1
  
  # 将 ER_factor 和 cluster 添加到数据框
  output_fc_top100 <- rbind(output_fc_top100, data.frame(fc_factor = fc_factor, cluster = cluster))
}

# 查看输出的数据框
print(output_fc_top100)

# 将数据框输出为 CSV 文件
write.csv(output_fc_top100, "fc_factor_top100_reedit.csv", row.names = FALSE)




output_fc <- data.frame(avg_log2FC = numeric(0), cluster = numeric(0))
# 循环遍历列表元素
for (i in 1:length(getsigpower_results)) {
  Sig_power <- getsigpower_results[[i]]$Sig_power
  
  Sig_power <- Sig_power[order(-Sig_power$avg_log2FC),]
  fc_factor <- head(rownames(Sig_power), 20)
  cluster <- i - 1
  
  # 将 ER_factor 和 cluster 添加到数据框
  output_fc <- rbind(output_fc, data.frame(fc_factor = fc_factor, cluster = cluster))
}

# 查看输出的数据框
print(output_fc)

# 将数据框输出为 CSV 文件
write.csv(output_fc, "fc_factor_reedit.csv", row.names = FALSE)



output_msi <- data.frame(msi_factor = numeric(0), cluster = numeric(0))
# 循环遍历列表元素
for (i in 1:length(getsigpower_results)) {
  Sig_power <- getsigpower_results[[i]]$Sig_power
  
  Sig_power <- Sig_power[order(-Sig_power$MSI),]
  msi_factor <- head(rownames(Sig_power), 20)
  cluster <- i - 1
  
  # 将 ER_factor 和 cluster 添加到数据框
  output_msi <- rbind(output_msi, data.frame(msi_factor = msi_factor, cluster = cluster))
}

# 查看输出的数据框
print(output_msi)

# 将数据框输出为 CSV 文件
write.csv(output_msi, "msi_factor_reedit.csv", row.names = FALSE)




# # 网络分析
gene_network_MI1 <- getsigpower_results[[1]]$G
gene_network_MI2 <- getsigpower_results[[2]]$G

# getsigpower_results[[6]]$ER_factor
# eg_mat <- getsigpower_results[[6]]$sparse_corr_mat

gene_network_MI <- gene_network_MI2

plot(gene_network_MI, vertex.label.color = "maroon",
     vertex.label.dist = 0.8,
     vertex.label.cex = 0.8,
     vertex.size = V(gene_network_MI)$size*2,
     edge.width = E(gene_network_MI)$width,
     edge.color = "grey69",
     layout = layout_with_kk(gene_network_MI),
     main = "Gene Co-exp Network",
     family = "Times New Roman")



# 创建一个空的网络列表 network_list
network_list <- vector("list", 2)
# 使用循环或者lapply函数将 getsigpower_results[[k]]$G 合并到 network_list 中
for (i in 1:2) {
  network_list[[i]] <- getsigpower_results[[i]]$G
  
  # 为每个节点添加序号属性
  V(network_list[[i]])$index <- i
}

# 初始化重叠节点为一个空列表
overlapping_nodes <- list()
# 逐对计算网络之间的重叠节点
for (i in 1:(length(network_list)-1)) {
  for (j in (i + 1):length(network_list)) {
    network1 <- network_list[[i]]
    network2 <- network_list[[j]]
    overlapping <- intersect(V(network1)$name, V(network2)$name)
    overlapping_nodes[[length(overlapping_nodes) + 1]] <- overlapping
  }
}

# 将list中所有的元素合并为统一的向量
commom_nodes <- unlist(overlapping_nodes)

# 创建一个空的多网络
merged_multi_network <- graph.empty()
# 合并网络列表中的所有网络
for (i in 1:length(network_list)) {
  current_network <- network_list[[i]]
  
  # 合并当前网络和已合并的多网络
  if (i == 1) {
    merged_multi_network <- current_network
  } else {
    merged_multi_network <- graph.union(merged_multi_network, current_network)
  }
}

# 计算平均节点大小
average_size <- rep(0, length(V(merged_multi_network)))
for (i in 1:length(network_list)) {
  current_network <- network_list[[i]]
  current_size <- V(current_network)$size
  
  average_size <- ifelse(V(merged_multi_network)$name %in% commom_nodes, current_size, current_size)
}

# 更新节点的颜色和大小属性
V(merged_multi_network)$color <- ifelse(V(merged_multi_network)$name %in% commom_nodes, "violet", "darkblue")
V(merged_multi_network)$size <- average_size

# 定义Force-Directed布局函数
layout <- layout_with_fr(merged_multi_network)
# 绘制合并后的多网络
plot(merged_multi_network,
     vertex.label.color = ifelse(V(merged_multi_network)$name %in% commom_nodes, "red", "black"),
     vertex.label.dist =1.25,
     vertex.label.cex = 0.5,
     vertex.size = V(merged_multi_network)$size*0.3,
     edge.width = E(merged_multi_network)$width,
     edge.color = "black",
     layout = layout,  # 使用Force-Directed布局
     main = "Merged Multi-Network",
     family = "Times New Roman")


# BiocManager::install("RCy3")
library(RCy3)
cytoscapePing()

createNetworkFromIgraph(merged_multi_network, "naiveESC Network")





##############  main-end 
getwd()
Sys.time()

.libPaths()


