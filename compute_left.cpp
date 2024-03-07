#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]

double compute_left(Eigen::VectorXd z, std::vector<Eigen::MatrixXd>& C, Eigen::SparseMatrix<double>& Sigma) {
  int N = z.size();
  double left_sum = 0;
  
  // Calcola il determinante di Sigma utilizzando la fattorizzazione LU
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.analyzePattern(Sigma);
  solver.factorize(Sigma);
  double det_Sigma = solver.determinant();

  std::vector<double> det_vec(N);
  for (int i = 0; i < N; ++i) {
    // Calcola il determinante di (0.5 * C[i])
    double det_C = (0.5 * C[i]).determinant();
    // std::cout << "Det C: " << det_C << std::endl;
    double log_detC = std::log(det_C);
    // std::cout << "LogDet C: " << log_detC << std::endl;
    det_vec[i] = log_detC - log(det_Sigma);
    
  }
  
  // std::cout << "Vector elements: ";
  // for (int i = 0; i < det_vec.size(); ++i) {
  //   std::cout << det_vec[i] << " ";
  // }
  // std::cout << std::endl;
  
  // Utilizza dot per il prodotto scalare
  left_sum = z.dot(Eigen::VectorXd::Map(det_vec.data(), det_vec.size()));
  
  return left_sum;
}
