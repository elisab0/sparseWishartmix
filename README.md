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

# Simulation Study README

This repository contains R code for a simulation study that evaluates the performance of the `compute_partition` function in estimating a partitioned covariance matrix.

## Objective

The primary goal of this simulation study is to assess the accuracy of the partitioning algorithm under various scenarios. The algorithm aims to partition a set of covariance matrices corresponding to different groups, considering partial observability of group memberships.

## Usage

To run the simulation, execute the provided R scripts:

1. **gen_sparse_posdef_mat_erdos_renyi.R**: Generates sparse positive definite matrices using an Erdos-Renyi model.
2. **compute_partition.R**: Implements the main function `compute_partition` and executes the simulation study.

Ensure that the required libraries (`Matrix`, `covglasso`, `stats`, `LaplacesDemon`, `igraph`, `mixggm`, `matrixStats`, `caret`, `mclust`, `ggplot2`, `e1071`, `Rcpp`) are installed before running the code.

## Simulation Parameters

Adjust the simulation parameters such as the number of simulations (`simulations`), dimensionality (`p`), the number of groups (`g`), and sample size (`N`) according to your experimental design.

## Results

The simulation produces results including F1 scores, Adjusted Rand Index (ARI), Mean Squared Error (MSE), and other metrics for each run. Results are saved in the `simulation_results_diff_sparse123_grid.rds` file.

## Analysis and Visualization

The R script includes code for visualizing F1 scores, MSE, ARI, and shrinkage values across simulations. The generated plots help analyze the algorithm's performance under different scenarios.


