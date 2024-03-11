#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

double compute_left(Eigen::VectorXd z, std::vector<Eigen::MatrixXd>& C, Eigen::SparseMatrix<double>& Sigma) {
  int N = z.size();
  double left_sum = 0;
  
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.analyzePattern(Sigma);
  solver.factorize(Sigma);
  double det_Sigma = solver.determinant();

  std::vector<double> det_vec(N);
  for (int i = 0; i < N; ++i) {
    double det_C = (0.5 * C[i]).determinant();
    double log_detC = std::log(det_C);
    det_vec[i] = log_detC - log(det_Sigma);
    
  }
  left_sum = z.dot(Eigen::VectorXd::Map(det_vec.data(), det_vec.size()));
  
  return left_sum;
}
