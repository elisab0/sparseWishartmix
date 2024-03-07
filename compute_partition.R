compute_partition <- function(C, g, p, Z_matrix, nu, shrinkage = .001, max_iter = 500, toll = 1e-6) {
  N <- length(C)  
  
  P <- rep(0, N)
  
  beta <- rep(0,g)
  
  Sigma_mat <- matrix(0,p,p)
  Sigma <- list()
  for(i in 1:g){
    Sigma[[i]] <- Sigma_mat
  }
  
  t <- 0
  
  obj_fun <- -Inf
  obj_fun_old <- -Inf
  
  BIC <- 0
  val = 0 
  
  criterion <- TRUE
  
  condition <- TRUE
  
  sourceCpp("equation_nu.cpp")
  sourceCpp("compute_left.cpp")
  
  while (criterion) {
    if (t != 0){
      for (i in 1:N){
        
        log_values <- numeric(g)
        
        for (j in 1:g) {
          log_values[j] <- log(beta[j]) + dwishart(Omega = C[[i]], nu = nu[j], S = Sigma[[j]], log = TRUE)
        }
        
        max_log_value <- max(log_values)
        exp_log_values <- exp(log_values - max_log_value)
        
        Z_matrix[i,] <- exp_log_values / sum(exp_log_values)
        
      }
    }
    
    # M-Step: Aggiorna n_k - risolvi numericamente, aggiorna Sigma_k, aggiorna p_k
    for (k in 1:g) {
      
      Z_sum = sum(Z_matrix[ ,k])
      
      beta[k] <- Z_sum / N
      
      C_new <- array(0, dim = c(p, p, N))
      
      for(i in 1:N){
        
        C_new[,,i] = Z_matrix[i,k] * C[[i]]
        
      }
      
      
      Sample_S <- apply(C_new,c(1,2),sum)/(Z_sum*nu[k])
      
      Sigma_sparse <- covglasso::covglasso(S = Sample_S,
                                           n = (Z_sum*nu[k]), # it should not matter
                                           lambda = rep(shrinkage, p))$sigma
      
      sigma_old <- matrix(0,p,p)
      
      criterion_sigma <- TRUE
      iter <- 0
      
      while (criterion_sigma) {
        
        iter <- iter + 1
        
        Sigma_sparse <- as(Sigma_sparse, "dgCMatrix")
        
        left_sum <- compute_left(z = Z_matrix[, k], C = C, Sigma = Sigma_sparse)
        
        Sigma_sparse <- as.matrix(Sigma_sparse)
        
        if(equation_nu(x = p, z = Z_matrix[, k], left_sum = left_sum, p = p) < 0){
          v_MLE <- p
          break
        }
        else{
          v_MLE <- pracma::bisect(equation_nu, a = p, b = 500, z = Z_matrix[, k], left_sum = left_sum, p = p, maxiter = 100)$root
        } 
        
        Sample_S <- apply(C_new, c(1, 2), sum) / (Z_sum * v_MLE)
        
        Sigma_sparse <- covglasso::covglasso(S = Sample_S,
                                             n = (Z_sum * v_MLE),
                                             lambda = rep(shrinkage, p))$sigma
        
        criterion_sigma <- (norm(sigma_old - Sigma_sparse, "F") > toll) & (iter < 250)
        sigma_old <- Sigma_sparse
        
      }
      
      nu[k] <- v_MLE
      
      Sigma[[k]] <- Sigma_sparse
      
    }
    
    # likelihood 
    dens <- matrix(NA, N, g)
    
    for (k in 1:g) {
      
      for (i in 1:N) {
        
        dens[i, k] <-
          
          dwishart(C[[i]],
                   
                   nu = nu[k],
                   
                   S = Sigma[[k]],
                   
                   log = TRUE)
        
      }
      
    }
    
    
    denspro <- sweep(dens, 2, log(beta), "+")
    
    zMax <- apply(denspro, 1, max)
    
    loghood <- zMax + log(rowSums(exp(denspro - zMax)))
    
    likelihood <- sum(loghood)
    
    penalty <- sum(sapply(1:g, function(k) shrinkage * sum(abs(Sigma[[k]]))))
    
    obj_fun <- likelihood - penalty
    
    # cat(obj_fun,"\n")
    
    criterion <- (abs(obj_fun - obj_fun_old) > toll) & (t < max_iter)
    
    obj_fun_old <- obj_fun
    
    print(t)
    
    t = t+1 
    
  }
  
  P <- apply(Z_matrix, 1, which.max)
  
  for(k in 1:g){
    val = val + p + sum(Sigma[[k]][-(1:nrow(Sigma[[k]])), ] != 0)/2
  }
  
  dof = g * (2+0.5*(val*(val+1))) - 1 
  
  BIC <- -2 * likelihood + dof * log(N)
  
  result = list(partition = P, Sigma = Sigma, nu = nu, BIC = BIC, likelihood = likelihood, penalty = penalty, shrinkage = shrinkage)
  
  return(result)
  
}
