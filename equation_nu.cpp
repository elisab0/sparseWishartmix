#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/math/special_functions/digamma.hpp>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double equation_nu(double x, Eigen::VectorXd z, double left_sum, int p) {
  int N = z.size();
  double psi = 0;
  
  for (int j = 1; j <= p; ++j) {
    psi += boost::math::digamma(0.5 * (x - j + 1));  // Utilizza la funzione di digamma di boost::math
  }
  
  double right_sum = z.dot(Eigen::VectorXd::Constant(N, psi));  // Utilizza dot per il prodotto scalare
  
  double result = left_sum - right_sum;
  
  return result;
}
