gen_sparse_posdef_mat_erdos_renyi <- function(prob_of_connections, V) {
  
  sample_graph <- igraph::as_adjacency_matrix(igraph::sample_gnp(n = V,
                                                                 
                                                                 p = prob_of_connections), sparse = FALSE)
  
  starting_matrix <- matrix(.9, nrow = V, ncol = V)
  
  diag(starting_matrix) <- 1
  
  mixggm:::icf(
    
    Sigma = diag(diag(starting_matrix), nrow = V, ncol = V),
    
    starting_matrix,
    
    sample_graph,
    
    10,
    
    mixggm::controlICF()$tol,
    
    mixggm::controlICF()$itMax,
    
    FALSE,
    
    FALSE,
    
    0
    
  )$sigma
  
}