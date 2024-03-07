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

source("C:\\Users\\ELISA\\Documents\\POLIMI\\TESI\\code\\gen_sparse_posdef_mat_erdos_renyi.R")
source("C:\\Users\\ELISA\\Documents\\POLIMI\\TESI\\code\\compute_partition.R")

set.seed(123)

simulations <- 50

simulation_results <- vector("list", length = simulations)

p = 20
g = 3 
N = 100

true_partition = rep(0,N*g)
true_partition[1:100] <-1
true_partition[101:200] <-2
true_partition[201:300] <-3

calculate_f1_score <- function(estimated_matrix, true_matrix) {
  # Elementi stimati come non zero e effettivamente non zero nella true matrix
  tp <- sum(estimated_matrix != 0 & true_matrix != 0)
  
  # Elementi stimati come zero ma sono non zero nella true matrix
  fp <- sum(estimated_matrix == 0 & true_matrix != 0)
  
  # Elementi stimati come non zero ma sono zero nella true matrix
  fn <- sum(estimated_matrix != 0 & true_matrix == 0)
  
  result <- tp/(tp+0.5*(fp+fn))
  
  return(result)
}


Sigma = list()

sparseness = c(0.2,0.4,0.6)
# sparseness = c(0.5,0.5,0.5)

for(it in 1:g){
  Sigma[[it]] <- gen_sparse_posdef_mat_erdos_renyi(sparseness[it], p)
}

nu <- c(20, 30, 40)


for (sim in 1:simulations){
  
  C1 <- rWishart(n = N, df = nu[1], Sigma = Sigma[[1]])
  
  C2 <- rWishart(n = N, df = nu[2], Sigma = Sigma[[2]])
  
  C3 <- rWishart(n = N, df = nu[3], Sigma = Sigma[[3]])
  
  arrayc1 <- array(C1,dim = c(p, p, N))
  
  C1l = vector("list", length = dim(arrayc1)[3])
  
  for (i in 1:dim(arrayc1)[3]) {
    C1l[[i]] <- arrayc1[,,i]
  }
  
  arrayc2 <- array(C2,dim = c(p, p, N))
  
  C2l = vector("list", length = dim(arrayc2)[3])
  
  for (i in 1:dim(arrayc2)[3]) {
    C2l[[i]] <- arrayc2[,,i]
  }
  
  arrayc3 <- array(C3,dim = c(p, p, N))
  
  C3l = vector("list", length = dim(arrayc3)[3])
  
  for (i in 1:dim(arrayc3)[3]) {
    C3l[[i]] <- arrayc3[,,i]
  }
  
  C <- c(C1l,C2l,C3l)
  
  Z_matrix <- matrix(0, length(C), g)
  
  Z_matrix[1:100,1] <-1
  Z_matrix[101:200,2] <-1
  Z_matrix[201:300,3] <-1
  
  vec_rand <- sample(1:300,180)
  
  for(i in vec_rand){
    if (i <= 100){
      k <- sample(c(2,3),1)
      Z_matrix[i,1] = 0
      Z_matrix[i,k] = 1
    }
    if (100 < i & i <= 200){
      k <- sample(c(1,3),1)
      Z_matrix[i,2] = 0
      Z_matrix[i,k] = 1
    }
    if (200 < i & i <= 300){
      k <- sample(c(1,2),1)
      Z_matrix[i,3] = 0
      Z_matrix[i,k] = 1
    }
  }
  
  shrinkage_values <- seq(0.001, 0.025, by = 0.006)
  
  results <- list()
  
  for (shrinkage in shrinkage_values) {
    
    partition <- compute_partition(C, g, p, Z_matrix, nu, shrinkage = shrinkage)
    
    results[[as.character(shrinkage)]] <- partition
  }
  
  best_index <- which.min(sapply(results, function(res) res$BIC))
  
  best_partition <- results[[best_index]]
  
  matching_clusters <-matchClasses(table(true_partition, best_partition$partition))

  # Per ogni simulazione
  
  simulation_results[[sim]] <- list(
    best_partition = best_partition,
    F1_score = c(calculate_f1_score(best_partition$Sigma[[matching_clusters[1]]], Sigma[[1]]),  
                    calculate_f1_score(best_partition$Sigma[[matching_clusters[2]]], Sigma[[2]]), 
                    calculate_f1_score(best_partition$Sigma[[matching_clusters[3]]], Sigma[[3]])),  
    ARI = adjustedRandIndex(best_partition$partition, true_partition), 
    MSE = c(norm(Sigma[[1]]-best_partition$Sigma[[matching_clusters[1]]], "F"), 
            norm(Sigma[[2]]-best_partition$Sigma[[matching_clusters[2]]], "F"),
            norm(Sigma[[3]]-best_partition$Sigma[[matching_clusters[3]]], "F"))  
  )
  
  print(sim)
  
}

# Salva i risultati della simulazione in un file binario
saveRDS(simulation_results, "simulation_results_diff_sparse123_grid.rds")

simulation_results = simulation_results_same_sparse123_grid

all_nu1 <- numeric(simulations)
all_nu2 <- numeric(simulations)
all_nu3 <- numeric(simulations)

for (i in 1:simulations){
  all_nu1[i] <- simulation_results[[i]]$best_partition$nu[1]
  all_nu2[i] <- simulation_results[[i]]$best_partition$nu[2]
  all_nu3[i] <- simulation_results[[i]]$best_partition$nu[3]
}

MSE_nu1 <- numeric(simulations)
MSE_nu2 <- numeric(simulations)
MSE_nu3 <- numeric(simulations)

# Calcola il MSE per i diversi cluster rispetto ai valori veri di nu
MSE_nu1 <- (all_nu1 - nu[1])^2
MSE_nu2 <- (all_nu2 - nu[2])^2
MSE_nu3 <- (all_nu3 - nu[3])^2

all_MSE_values1 <- numeric(simulations)
all_MSE_values2 <- numeric(simulations)
all_MSE_values3 <- numeric(simulations)

all_ARI_scores <- numeric(simulations)

all_F1_scores1 <- numeric(simulations)
all_F1_scores2 <- numeric(simulations)
all_F1_scores3 <- numeric(simulations)

for (i in 1:simulations) {
  all_MSE_values1[i] <- simulation_results[[i]]$MSE[1]
  all_MSE_values2[i] <- simulation_results[[i]]$MSE[2]
  all_MSE_values3[i] <- simulation_results[[i]]$MSE[3]
  
  all_ARI_scores[i] <- simulation_results[[i]]$ARI
  
  all_F1_scores1[i] <- simulation_results[[i]]$F1_score[1]
  all_F1_scores2[i] <- simulation_results[[i]]$F1_score[2]
  all_F1_scores3[i] <- simulation_results[[i]]$F1_score[3]
}

# means_MSE <- c(mean(all_MSE_values1), mean(all_MSE_values2), mean(all_MSE_values3))
# means_ARI = mean(all_ARI_scores) 
# means_F1score <- c(mean(all_F1_scores1), mean(all_F1_scores2), mean(all_F1_scores3))


#### PLOTS #####################################################################
results_df <- data.frame(
  MSE_nu1,
  MSE_nu2,
  MSE_nu3,
  all_F1_scores1,
  all_F1_scores2,
  all_F1_scores3,
  all_MSE_values1,
  all_MSE_values2,
  all_MSE_values3,
  all_ARI_scores
)

# Grafico boxplot per F1_score
ggplot(results_df, aes(x = all_F1_scores1, y = factor(1))) +
  geom_boxplot() +
  labs(title = "F1 Score - Cluster 1")

ggplot(results_df, aes(x = all_F1_scores2, y = factor(1))) +
  geom_boxplot() +
  labs(title = "F1 Score - Cluster 2")

ggplot(results_df, aes(x = all_F1_scores3, y = factor(1))) +
  geom_boxplot() +
  labs(title = "F1 Score - Cluster 3")

# Grafico boxplot per MSE
ggplot(results_df, aes(x = all_MSE_values1, y = factor(1))) +
  geom_boxplot() +
  labs(title = "MSE - Cluster 1")

ggplot(results_df, aes(x = all_MSE_values2, y = factor(2))) +
  geom_boxplot() +
  labs(title = "MSE - Cluster 2")

ggplot(results_df, aes(x = all_MSE_values3, y = factor(3))) +
  geom_boxplot() +
  labs(title = "MSE - Cluster 3")


# Grafico boxplot per MSE Nu
ggplot(results_df, aes(x = MSE_nu1, y = factor(1))) +
  geom_boxplot() +
  labs(title = "MSE Nu - Cluster 1")

ggplot(results_df, aes(x = MSE_nu2, y = factor(2))) +
  geom_boxplot() +
  labs(title = "MSE Nu - Cluster 2")

ggplot(results_df, aes(x = MSE_nu3, y = factor(3))) +
  geom_boxplot() +
  labs(title = "MSE Nu - Cluster 3")

# Grafico boxplot per ARI
ggplot(results_df, aes(x = all_ARI_scores, y = factor(1))) +
  geom_boxplot() +
  labs(title = "Adjusted Rand Index")

# # Grafico boxplot per nu
# ggplot(results_df, aes(x = all_nu1, y = factor(1))) +
#   geom_boxplot() +
#   labs(title = "Nu - Cluster 1")
# 
# ggplot(results_df, aes(x = all_nu2, y = factor(1))) +
#   geom_boxplot() +
#   labs(title = "Nu - Cluster 2")
# 
# ggplot(results_df, aes(x = all_nu3, y = factor(1))) +
#   geom_boxplot() +
#   labs(title = "Nu - Cluster 3")
# 
# 
# 
# plot(all_F1_scores1, type = "p", col = "blue", lwd = 2, xlab = "Simulations", ylab = "F1 Score", main = "Plot F1 Scores-Cluster 1", ylim = c(0, 1))
# abline(h = means_F1score[1], col = "red", lty = 1)
# legend("topright", legend = c("F1 Scores", "Mean F1"), col = c("blue", "red"), lty = c(1, 2), lwd = c(2, 1), cex = 0.8)
# 
# plot(all_F1_scores2, type = "p", col = "blue", lwd = 2, xlab = "Simulations", ylab = "F1 Score", main = "Plot F1 Scores-Cluster 2", ylim = c(0, 1))
# abline(h = means_F1score[2], col = "red", lty = 1)
# legend("topright", legend = c("F1 Scores", "Mean F1"), col = c("blue", "red"), lty = c(1, 2), lwd = c(2, 1), cex = 0.8)
# 
# plot(all_F1_scores3, type = "p", col = "blue", lwd = 2, xlab = "Simulations", ylab = "F1 Score", main = "Plot F1 Scores-Cluster 3", ylim = c(0, 1))
# abline(h = means_F1score[3], col = "red", lty = 1)
# legend("topright", legend = c("F1 Scores", "Mean F1"), col = c("blue", "red"), lty = c(1, 2), lwd = c(2, 1), cex = 0.8)
# 
# plot(all_MSE_values1, type = "p", col = "blue", lwd = 2, xlab = "Simulations", ylab = "MSE", main = "Plot MSE-Cluster 1", ylim = c(0, 10))
# abline(h = means_MSE[1], col = "red", lty = 1)
# legend("topright", legend = c("MSE", "Mean MSE"), col = c("blue", "red"), lty = c(1, 2), lwd = c(2, 1), cex = 0.8)
# 
# plot(all_MSE_values2, type = "p", col = "blue", lwd = 2, xlab = "Simulations", ylab = "MSE", main = "Plot MSE-Cluster 2", ylim = c(0, 10))
# abline(h = means_MSE[2], col = "red", lty = 1)
# legend("topright", legend = c("MSE", "Mean MSE"), col = c("blue", "red"), lty = c(1, 2), lwd = c(2, 1), cex = 0.8)
# 
# plot(all_MSE_values3, type = "p", col = "blue", lwd = 2, xlab = "Simulations", ylab = "MSE", main = "Plot MSE-Cluster 3", ylim = c(0, 10))
# abline(h = means_MSE[3], col = "red", lty = 1)
# legend("topright", legend = c("MSE", "Mean MSE"), col = c("blue", "red"), lty = c(1, 2), lwd = c(2, 1), cex = 0.8)
# 
# plot(all_ARI_scores, type = "o", col = "blue", lwd = 2, xlab = "Simulations", ylab = "ARI Scores", main = "Plot ARI Scores")
# abline(h = means_ARI, col = "red", lty = 1)
# legend("topright", legend = c("ARI Scores", "Mean ARI"), col = c("blue", "red"), lty = c(1, 2), lwd = c(2, 1), cex = 0.8)
# 
# plot(all_nu1, type = "p", col = "blue", lwd = 2, xlab = "Simulations", ylab = "Nu", main = "Plot Nu-Cluster 1")
# abline(h = nu1_means, col = "red", lty = 1)
# legend("topright", legend = c("Nu", "Mean Nu"), col = c("blue", "red"), lty = c(1, 2), lwd = c(2, 1), cex = 0.8)
# 
# plot(all_nu2, type = "p", col = "blue", lwd = 2, xlab = "Simulations", ylab = "Nu", main = "Plot Nu-Cluster 2")
# abline(h = nu2_means, col = "red", lty = 1)
# legend("topright", legend = c("Nu", "Mean Nu"), col = c("blue", "red"), lty = c(1, 2), lwd = c(2, 1), cex = 0.8)
# 
# plot(all_nu3, type = "p", col = "blue", lwd = 2, xlab = "Simulations", ylab = "Nu", main = "Plot Nu-Cluster 3")
# abline(h = nu3_means, col = "red", lty = 1)
# legend("topright", legend = c("Nu", "Mean Nu"), col = c("blue", "red"), lty = c(1, 2), lwd = c(2, 1), cex = 0.8)
# 

# Estrai i valori di shrinkage da ciascuna best_partition
shrinkage_values_used <- unlist(lapply(simulation_results, function(sim) sim$best_partition$shrinkage))

# Conta la frequenza di ciascun valore di shrinkage
shrinkage_freq <- table(shrinkage_values_used)

# Crea un dataframe con i risultati
shrinkage_df <- data.frame(shrinkage = as.numeric(names(shrinkage_freq)),
                           frequency = as.numeric(shrinkage_freq))

# Visualizza un istogramma
hist_plot <- ggplot(shrinkage_df, aes(x = shrinkage, y = frequency)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  labs(title = "Frequenza dei Valori di Shrinkage", x = "Shrinkage", y = "Frequenza") +
  theme_minimal()

print(hist_plot)

# Visualizza un grafico a torta
pie_chart <- ggplot(shrinkage_df, aes(x = "", y = frequency, fill = factor(shrinkage))) +
  geom_bar(stat = "identity", color = "white", width = 1) +
  coord_polar("y") +
  labs(title = "Shrinkage Values Distribution", x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank())

print(pie_chart)




