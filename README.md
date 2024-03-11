# compute_partition
The `compute_partition` function is designed to perform partitioning using a specified set of parameters. This function employs an E-M iterative optimization process to estimate parameters such as weights, covariance matrices, and degrees of freedom. The primary goal is to maximize the likelihood of the observed data while penalizing the complexity of the model.

## Function Parameters

- `C`: A list of covariance matrices.
- `g`: The number of groups.
- `p`: The dimension of the covariance matrices.
- `Z_matrix`: A matrix to store intermediate results during the optimization process.
- `nu`: A vector representing the degrees of freedom for each group.
- `shrinkage`: A parameter controlling the amount of regularization applied to the covariance matrices (default is 0.001).
- `max_iter`: The maximum number of iterations for the optimization process (default is 500).
- `toll`: Convergence tolerance for the iterative optimization process (default is 1e-6).

## Iterative Optimization Process

The function utilizes an iterative algorithm to update parameters until convergence or reaching the maximum number of iterations. The steps include updating weights, covariance matrices, and degrees of freedom for each group.

### External Dependencies

- `equation_nu.cpp`: C++ code containing an equation related to degrees of freedom.
- `compute_left.cpp`: C++ code providing a computation related to degrees of freedom, to give in input to equation_nu.cpp.
  The choice of writing it in C++ is due to computational speed. 

## Outputs

The function returns a list containing the following results:

- `partition`: A vector indicating the assignment of each observation to a specific group.
- `Sigma`: A list of estimated covariance matrices for each group.
- `nu`: A vector of estimated degrees of freedom for each group.
- `BIC`: The Bayesian Information Criterion calculated based on the likelihood and model complexity.
- `likelihood`: The log-likelihood of the observed data given the model.
- `penalty`: The penalty term accounting for the complexity of the model.
- `shrinkage`: The specified shrinkage parameter used in the optimization process.

