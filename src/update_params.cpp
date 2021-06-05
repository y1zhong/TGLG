#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]


arma::mat update_test(arma::vec a){
  mat v = randu<mat>(3,3);
  mat b =randu<mat>(3,3);
  mat c =  v * b*a;
  return(c);
}
arma::vec dappro2(arma::vec x, double epison = 10^-8) {
  vec temp = 2/(datum::pi*(epison+pow(x,2)/epison));
  return(temp);
}

// [[Rcpp::export]]
List update_gamma(arma::vec gamma, arma::uvec alpha, arma::vec y,
                        arma::mat X, arma::mat eigmat,
                       double lambda, double sigmae,
                       double sigmagamma, double taugamma){
  vec currgamma = gamma;
  vec currhgamma = pow(currgamma,2)-pow(lambda,2);
  uvec currbeta = alpha % (abs(currgamma)>lambda);
  vec yerr = y - X*currbeta;
  vec currdiff = X.t()*yerr;
  vec currdgamma = dappro2(currhgamma) % currgamma;
  vec eigmat_currgamma = eigmat*currgamma;
  vec currpost = sum(pow(yerr,2))/sigmae*0.5+0.5*currgamma.t()*(eigmat * currgamma)/sigmagamma;
  vec newgammamean = currgamma + taugamma/2*(currdiff*gamma.t()*currdgamma/sigmae - eigmat*currgamma/sigmagamma);
  vec newrnorm = rnorm(gamma.n_elem,0,sqrt(taugamma));
  vec newgamma = newgammamean+newrnorm;
  vec newhgamma = pow(newgamma,2)-pow(lambda,2);
  uvec newbeta=alpha%(abs(newgamma)>lambda);
  vec newyerr = y - X*newbeta;
  vec newdiff = X.t()*newyerr;
  //Rcout << "halfway " << std::endl;
  vec newdgamma = dappro2(newhgamma)%newgamma;
  vec currgammamean = newgamma + taugamma/2*(newdiff*gamma.t()*newdgamma/sigmae - eigmat*newgamma/sigmagamma);
  vec newpost=sum(pow(newyerr,2))/sigmae*0.5+0.5*newgamma.t()*(eigmat*newgamma)/sigmagamma;
  vec curradiff = newgamma - newgammamean;
  vec currden = 0.5*sum(pow(curradiff,2))/taugamma;
  vec newadiff = currgamma - currgammamean;
  vec newden = 0.5* sum(pow(newadiff,2))/taugamma;
  vec postdiff = currpost + currden - newpost - newden;
  List L =List::create(Named("postdiff")=postdiff,Named("newgamma")= newgamma);
  return(L);
}


arma::vec update_alpha(vec curr.gamma, vec lambda, arma::mat X) {
  premat = XA)/sigmae;
  diag(premat) <- diag(premat) + (1/sigma.alpha);
  premat_chol <- chol(premat);
  mu <- crossprod(XA, y)/sigmae;
  b <- backsolve(premat_chol, mu, transpose=T);
  Z <- rnorm(actset_len, 0, 1);
  alpha[actset] <- backsolve(premat_chol, Z+b);
  //return(alpha);
}

//arma::vec update_sigma(vec curr.gamma, vec lambda, arma::mat X) {
  //return(sigma);
//}

//arma::vec update_epsilon(vec curr.gamma, vec lambda, arma::mat X) {
  //return(epsilon);
//}

//arma::vec update_lambda(vec curr.gamma, vec lambda, arma::mat X) {
  //return(sigma);
//}