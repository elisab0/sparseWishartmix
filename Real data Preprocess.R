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
library(gplots)
library(grid)
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
x11(width = 8, height = 7)

heatmap.2(as.matrix(D_average),
          col = custom_palette(length(color_scale)),
          key = FALSE,
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          trace = "none",
          main = "Heatmap White Fibers",
          labRow = FALSE,  
          labCol = FALSE)


legend("left", legend = pretty(color_scale), fill = custom_palette(length(color_scale)), title = "Valori")

quantile(D_average, probs = seq(from = 0, to = 1, by = 0.05), na.rm = T) 

# Creare una matrice binarizzata
D_binary <- ifelse(D_average < 0.1, 0, 1)

x11(width = 8, height = 7)


heatmap.2(as.matrix(D_binary),
          col = c("lightblue", "darkred"),
          key = FALSE,
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          trace = "none",
          main = "Heatmap Binary White Fibers",
          labRow = FALSE, 
          labCol = FALSE)


legend("left", legend = c("No fibers", "Fibers"), fill = c("lightblue", "darkred"), title = "Connections")


connections_count <- rowSums(D_binary, na.rm = TRUE)

mean(connections_count)
quantile(connections_count, probs = seq(from = 0, to = 1, by = 0.05), na.rm = T) 

ind_zero = which(connections_count <= 30)

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
heatmap.2(as.matrix(W_average),
          col = custom_palette(length(color_scale)),
          key = FALSE,
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          trace = "none",
          main = "Heatmap mean fMRI",
          labRow = FALSE,  # Imposta a FALSE per rimuovere le etichette delle righe
          labCol = FALSE)


x11()
plot.new()
legend("center", legend = legend_labels, fill = color_scale, title = "Valori", xpd = TRUE)

# Matrice dei range
limiti <- c(0.155, 0.324, 0.493, 0.662, 0.831, 1)

W_range <- matrix(0, nrow = nrow(W_average), ncol = ncol(W_average))

W_range[W_average >= limiti[1] & W_average < limiti[2]] <- 1
W_range[W_average >= limiti[2] & W_average < limiti[3]] <- 2
W_range[W_average >= limiti[3] & W_average < limiti[4]] <- 3
W_range[W_average >= limiti[4] & W_average < limiti[5]] <- 4
W_range[W_average >= limiti[5] & W_average <= limiti[6]] <- 5

custom_palette <- colorRampPalette(c("lightblue", "lightgreen","pink", "orange", "darkred"))

x11(width = 8, height = 7)

heatmap.2(as.matrix(W_range),
          col = custom_palette(5),
          key = FALSE,
          dendrogram = "none",
          Rowv = FALSE,
          Colv = FALSE,
          trace = "none",
          main = "Heatmap Range Matrix",
          labRow = FALSE,  # Imposta a FALSE per rimuovere le etichette delle righe
          labCol = FALSE)

legend("left", legend = seq(1, 5), fill = custom_palette(5), title = "Valori")


cont1 <- rowSums(W_range == 1)
cont2 <- rowSums(W_range == 2)
cont3 <- rowSums(W_range == 3)
cont4 <- rowSums(W_range == 4)
cont5 <- rowSums(W_range == 5)

comb <- cbind(cont1, cont2, cont3, cont4, cont5)
