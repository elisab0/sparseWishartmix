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


D_list <- vector("list", 24)

for (i in 1:24) {
  mat <- D[, , i, 1]
  D_list[[i]] <- mat
}


empty_matrix_positions <- which(sapply(D_list, function(mat) sum(!is.na(mat)) == 0))

D_list <- D_list[-empty_matrix_positions]

D_average <- Reduce(`+`, D_list) / length(D_list)

# Unknown 
D_average <- D_average[-c(1,36), -c(1,36)]

rownames(D_average) <- NULL
colnames(D_average) <- NULL

color_scale <- seq(min(D_average[!is.na(D_average)]), 
                   max(D_average[!is.na(D_average)]), 
                   length.out = 11)

custom_palette <- colorRampPalette(c("lightblue", "blue", "white", "pink", "darkred"))

x11()
image(as.matrix(D_average), col = custom_palette(length(color_scale)), axes = FALSE, main = "Heatmap White Fibers")

legend("topright", legend = pretty(color_scale), fill = custom_palette(length(color_scale)), title = "Valori")

quantile(D_average, probs = seq(from = 0, to = 1, by = 0.05), na.rm = T) 

# Creare una matrice binarizzata
D_binary <- ifelse(D_average < 0.1, 0, 1)


x11()
image(as.matrix(D_binary), col = c("lightblue", "darkred"), axes = FALSE, main = "Heatmap Binary White Fibers")
legend("topright", legend = c("No fibers", "Fibers"), fill = c("lightblue", "darkred"), title = "Connections")

# Contare il numero di connessioni per ogni riga
connections_count <- rowSums(D_binary, na.rm = TRUE)

mean(connections_count)
quantile(connections_count, probs = seq(from = 0, to = 1, by = 0.05), na.rm = T) 

ind_zero = which(connections_count <= 30)

# Escludo le righe e colonne che hanno connection_counts == 0 

D_average = D_average[-ind_zero, -ind_zero]

################################################################################
# Tengo conto solo del primo scan perchè il secondo è presente solo per pochissimi soggetti  
W_list <- vector("list", 24)

for (i in 1:24) {
  
  mat <- W[, , i, 1]
  W_list[[i]] <- mat
  
}

empty_matrix_positions <- which(sapply(W_list, function(mat) sum(!is.na(mat)) == 0))

# Check and remove empty matrices from the list
W_list <- W_list[-empty_matrix_positions]

W_average <- Reduce(`+`, W_list) / length(W_list)

# Unknown 
W_average <- W_average[-c(1,36), -c(1,36)]

rownames(W_average) <- NULL
colnames(W_average) <- NULL

# Escludo le righe e colonne che hanno connection_counts == 0 

W_average = W_average[-ind_zero, -ind_zero]

custom_palette <- colorRampPalette(c("lightblue", "blue", "white", "pink", "darkred"))

num_colors <- 11
color_scale <- custom_palette(num_colors)

legend_labels <- seq(min(W_average, na.rm = TRUE), max(W_average, na.rm = TRUE), length.out = num_colors)

x11()
image(as.matrix(W_average), col = color_scale, axes = FALSE, main = "Heatmap Mean fMRI")

legend("topright", legend = legend_labels, fill = color_scale, title = "Valori")

# Matrice dei range
limiti <- c(0.155, 0.324, 0.493, 0.662, 0.831, 1)

W_range <- matrix(0, nrow = nrow(W_average), ncol = ncol(W_average))

W_range[W_average >= limiti[1] & W_average < limiti[2]] <- 1
W_range[W_average >= limiti[2] & W_average < limiti[3]] <- 2
W_range[W_average >= limiti[3] & W_average < limiti[4]] <- 3
W_range[W_average >= limiti[4] & W_average < limiti[5]] <- 4
W_range[W_average >= limiti[5] & W_average <= limiti[6]] <- 5

custom_palette <- colorRampPalette(c("lightblue", "lightgreen","pink", "orange", "darkred"))

x11()
image(as.matrix(W_range), col = custom_palette(5), axes = FALSE, main = "Heatmap Range Matrix")

legend("topright", legend = seq(1, 5), fill = custom_palette(5), title = "Valori")

cont1 <- rowSums(W_range == 1)
cont2 <- rowSums(W_range == 2)
cont3 <- rowSums(W_range == 3)
cont4 <- rowSums(W_range == 4)
cont5 <- rowSums(W_range == 5)

comb <- cbind(cont1, cont2, cont3, cont4, cont5)

# elim <- c()
# j = 1
# 
# for(i in 1:dim(comb)[1]){
#   if(comb[i,1] + comb[i,2] >= 50){
#     elim[j] = i
#     j = j+1
#   }
# }
# 
# elim = unique(elim)
# 
# W_average = W_average[-elim, -elim]

################################################################################
for(i in 1:length(W_list)){
  W_list[[i]] <- W_list[[i]][-c(1,36), -c(1,36)]
}


for(i in 1:length(W_list)){
  W_list[[i]] <- W_list[[i]][-ind_zero,-ind_zero]
}


# for(i in 1:length(W_list)){
#   W_list[[i]] <- W_list[[i]][-elim,-elim]
# }

p = 59
nu = c(60,60,60)
# shrinkage = 0.06

data <- sapply(W_list, function(mat) as.vector(mat))

data <- t(data)

b <- NULL
w <- NULL
for(k in 1:10){
  
  result.k <- kmeans(data, k)
  w <- c(w, sum(result.k$wit))
  b <- c(b, result.k$bet)
  
}

x11()
matplot(1:10, w/(w+b), pch='', xlab='clusters', ylab='within/tot', main='Choice of k', ylim=c(0,1))
lines(1:10, w/(w+b), type='b', lwd=2)

result.k <- kmeans(data, 2)

cluster_assegnati <- result.k$cluster

g = 2

Z_matrix <- matrix(0, length(W_list), g)

for (i in 1:length(cluster_assegnati)){
  Z_matrix[i, cluster_assegnati[i]] <- 1
}

non_posdef <- c()
j = 1

for (i in 1:length(W_list)){
  if(!is.positive.semidefinite(W_list[[i]])){
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

#shrinkage_values <- seq(10, 15, by = 1)

g_values <- seq(2, 10, by = 1)

results <- list()

# for (shrinkage in shrinkage_values) {
# 
#   partition <- compute_partition(C = W_list, g = g, p = p, Z_matrix = Z_matrix, nu = nu, shrinkage = shrinkage)
# 
#   results[[as.character(shrinkage)]] <- partition
# }

for (g in g_values) {

    partition <- compute_partition(C = W_list, g = g, p = p, Z_matrix = Z_matrix, nu = nu, shrinkage = 15)

    results[[as.character(g)]] <- partition
  }
  


best_index <- which.min(sapply(results, function(res) res$BIC))

best_partition <- results[[best_index]]

saveRDS(best_partition, "best_partition_70_10-15_g2_onlyind0.rds")
saveRDS(results, "results_70_10-15_g2_onlyind0.rds")

sum(best_partition$Sigma[[1]] != 0)
sum(best_partition$Sigma[[2]] != 0)
sum(best_partition$Sigma[[3]] != 0)

best_partition$nu[1]
best_partition$nu[2]
best_partition$nu[3]

sum(results$'10'$Sigma[[1]] != 0)
sum(results$'10'$Sigma[[2]] != 0)
sum(results$'10'$Sigma[[3]] != 0)

sum(results$'11'$Sigma[[1]] != 0)
sum(results$'11'$Sigma[[2]] != 0)
sum(results$'11'$Sigma[[3]] != 0)

sum(results$'12'$Sigma[[1]] != 0)
sum(results$'12'$Sigma[[2]] != 0)
sum(results$'12'$Sigma[[3]] != 0)

sum(results$'13'$Sigma[[1]] != 0)
sum(results$'13'$Sigma[[2]] != 0)
sum(results$'13'$Sigma[[3]] != 0)

sum(results$'14'$Sigma[[1]] != 0)
sum(results$'14'$Sigma[[2]] != 0)
sum(results$'14'$Sigma[[3]] != 0)

sum(results$'15'$Sigma[[1]] != 0)
sum(results$'15'$Sigma[[2]] != 0)
sum(results$'15'$Sigma[[3]] != 0)

# vedere se c'è differenza tra P e Z

table(best_partition$partition,cluster_assegnati)

x11()
plot(data, col = best_partition$partition)
