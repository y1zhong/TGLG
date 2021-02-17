#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::vec eigen_val(arma::mat X) {
  //int k = X.n_rows -1 ;
  vec eigval=eig_sym(X);
  return(eigval);
}