library(Matrix)
library(covglasso)
library(stats)
library(LaplacesDemon)
library(igraph)
library(mixggm)
library(matrixStats)
library(caret)
library(mclust)
library(ggplot2)
library(e1071)
library(Rcpp)
library(glasso)
library(corpcor)
source("C:\\Users\\ELISA\\Documents\\POLIMI\\TESI\\code\\compute_partition.R")

# Tengo conto solo del primo scan perchè il secondo è presente solo per pochissimi soggetti  
W_list <- vector("list", 24)

for (i in 1:24) {
  
  mat <- W[, , i, 1]
  W_list[[i]] <- mat
  
}

empty_matrix_positions <- which(sapply(W_list, function(mat) sum(!is.na(mat)) == 0))

W_list <- W_list[-empty_matrix_positions]

na_check <- sapply(W_list, function(mat) any(is.na(mat)))

p = 70

non_posdef <- c()
j = 1

for (i in 1:length(W_list)){
  if(!is.positive.definite(W_list[[i]])){
    non_posdef[j] <- i
    j = j+1
  }
}

for(i in 1:length(W_list)){
  if(!is.positive.definite(W_list[[i]])){
    W_list[[i]] <- make.positive.definite(W_list[[i]])
  }
}

non_posdef <- c()
j = 1

for (i in 1:length(W_list)){
  if(!is.positive.definite(W_list[[i]])){
    non_posdef[j] <- i
    j = j+1
  }
}

for ( i in 1:length(W_list)){
  if(!is.symmetric.matrix(W_list[[i]])){
    print(i)
    W_list[[i]] <- (W_list[[i]] + t(W_list[[i]]))/2
    print(is.symmetric.matrix(W_list[[i]]))
  }
}

################################################################################

g_values <- seq(2, 10, by = 1)

results <- list()

N = length(W_list)

W_tensor <- array(unlist(W_list), dim = c(p, p, N))
W_matrix <- matrix(W_tensor, nrow = p * p, ncol = N)

for (g in g_values) {
  
  cat("G:", g, "\n")

  result_clus <- cmeans(t(W_matrix), centers = g)
  
  Z_matrix <- result_clus$membership
  
  nu = rep(80,g)
  
  partition <- compute_partition(C = W_list, g = g, p = p, Z_matrix = Z_matrix, nu = nu)
  
  results[[as.character(g)]] <- partition
}

best_index <- which.min(sapply(results, function(res) res$BIC))

best_partition <- results[[best_index]]

saveRDS(results, "results_fin_FCM")

###############################################################################

shrinkage_values <- seq(0.05, 1.05, by = 0.2)

results_shrink <- list()

result_clus <- cmeans(t(W_matrix), centers = 2)

Z_matrix <- result_clus$membership

nu = rep(80,2)

for (shrinkage in shrinkage_values) {

   partition <- compute_partition(C = W_list, g = 2, p = p, Z_matrix = Z_matrix, nu = nu, shrinkage = shrinkage)

   results_shrink[[as.character(shrinkage)]] <- partition
}

best_index <- which.min(sapply(results_shrink, function(res) res$BIC))

best_partition <- results_shrink[[best_index]]

sum(best_partition$Sigma[[1]] != 0)
sum(best_partition$Sigma[[2]] != 0)

best_partition$nu[1]
best_partition$nu[2]

best_partition$partition

###############################################################################
nu = c(80,80)

g = 2 

result_clus <- cmeans(t(W_matrix), centers = g)

centroids <- t(result_clus$centers)
membership_probs <- result_clus$membership

Z_matrix <- membership_probs

partition <- compute_partition(C = W_list, g = g, p = p, Z_matrix = Z_matrix, nu = nu, shrinkage = 15)

S_binary1 <- ifelse(partition$Sigma[[1]] == 0, 0, 1)
S_binary2 <- ifelse(partition$Sigma[[2]] == 0, 0, 1)

x11()
image(as.matrix(S_binary1), col = c("white", "black"), axes = FALSE, main = "Sigma 1")

x11()
image(as.matrix(S_binary2), col = c("white", "black"), axes = FALSE, main = "Sigma 2")

sum(partition$Sigma[[1]] == 0) 
sum(partition$Sigma[[2]] == 0)

