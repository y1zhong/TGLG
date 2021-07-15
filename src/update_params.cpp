// #include <RcppArmadillo.h>
// using namespace Rcpp;
// using namespace arma;
// 
// #define MATH_PI        3.141592653589793238462643383279502884197169399375105820974
// Environment pkg = Environment::namespace_env("glmnet");
// Function f1 = pkg["cv.glmnet"];
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// arma::vec dappro2(arma::vec& x, double epison = pow(10,-8)) {
//   vec temp = 2/((MATH_PI)*(epison+square(x)/epison));
//   return(temp);
// }
// 
// List TGLG_continuous_rcpp(arma::mat X, vec y, arma::mat L=NULL, arma::vec beta_ini = NULL,
//                           int nsim=30000, int ntune=10000, int freqTune=100, 
//                           int nest=10000, int nthin=10, 
//                           double a_alpha=0.01,double b_alpha=0.01,
//                           double a_gamma=0.01,double b_gamma=0.01,
//                           double ae0=0.01,double be0=0.01,
//                           double lambda_lower=0, double lambda_upper=10,
//                           double emu=-5, double esd=3, 
//                           double prob_select=0.95){
//   int p = X.n_cols;
//   
//   if(is_na(L)){
//     L = mat(p, p);
//     
//   }
// 
//   double tau_gamma=0.01/p;
//   double tau_lambda=0.1;
//   double tau_epsilon = 0.5;
//   
//   //count number of acceptance during MCMC updates
//   int accept_gamma=0;
//   int accept_epsilon=0;
//   int accept_lambda =0;
//   
//   vec beta_elas=as<vec>(f1(X,y,Named("alpha", 1),Named("intercept", "FALSE")));
//   
// }
// 
// // [[Rcpp::depends(RcppArmadillo)]]
// // [[Rcpp::export]]
// List update_gamma(arma::vec gamma, arma::vec alpha, arma::vec y,
//                         arma::mat X, arma::mat eigmat,
//                        double lambda, double sigmae,
//                        double sigmagamma, double taugamma){
//   vec currgamma = gamma;
//   vec currhgamma = pow(currgamma,2)-pow(lambda,2);
//   uvec gammaidx = find(abs(currgamma)<=lambda);
//   vec currbeta = alpha;
//   currbeta.elem(gammaidx) = zeros(gammaidx.n_elem);
//   //Rcout << "currbeta " <<currbeta<< std::endl;
//   vec yerr = y - X*currbeta;
//   vec currdiff = X.t()*yerr;
//   vec currdgamma = dappro2(currhgamma) % currgamma;
//   vec eigmat_currgamma = eigmat*currgamma;
//   vec currpost = sum(pow(yerr,2))/sigmae*0.5+0.5*currgamma.t()*(eigmat * currgamma)/sigmagamma;
//   vec newgammamean = currgamma + taugamma/2*(currdiff*gamma.t()*currdgamma/sigmae - eigmat*currgamma/sigmagamma);
//   vec newrnorm = rnorm(gamma.n_elem,0,sqrt(taugamma));
//   vec newgamma = newgammamean+newrnorm;
//   vec newhgamma = pow(newgamma,2)-pow(lambda,2);
//   //uvec newbeta=alpha%(abs(newgamma)>lambda);
//   uvec newgammaidx = find(abs(newgamma)<=lambda);
//   vec newbeta = alpha;
//   newbeta.elem(newgammaidx) = zeros(newgammaidx.n_elem);
//   vec newyerr = y - X*newbeta;
//   vec newdiff = X.t()*newyerr;
//   //Rcout << "halfway " << std::endl;
//   vec newdgamma = dappro2(newhgamma)%newgamma;
//   vec currgammamean = newgamma + taugamma/2*(newdiff*gamma.t()*newdgamma/sigmae - eigmat*newgamma/sigmagamma);
//   vec newpost=sum(pow(newyerr,2))/sigmae*0.5+0.5*newgamma.t()*(eigmat*newgamma)/sigmagamma;
//   vec curradiff = newgamma - newgammamean;
//   vec currden = 0.5*sum(pow(curradiff,2))/taugamma;
//   vec newadiff = currgamma - currgammamean;
//   vec newden = 0.5* sum(pow(newadiff,2))/taugamma;
//   vec postdiff = currpost + currden - newpost - newden;
//   List L =List::create(Named("postdiff")=postdiff,Named("newgamma")= newgamma);
//   return(L);
// }
// 
// // [[Rcpp::export]]
// arma::vec update_alpha(arma::vec gamma, arma::vec alpha, arma::vec y,
//                        arma::mat X,
//                        double lambda, double sigmaalpha, double sigmae) {
//   uvec actset = find(abs(gamma)>lambda);
//   uvec nonactset = find(abs(gamma) <= lambda);
//   int actset_len = actset.n_elem;
//   int nonactset_len = nonactset.n_elem;
//   Function f("backsolve");
//   vec norm;
//   //Rcout << "element"<<nonactset;
//   //Rcout << "rnorm"<<rnorm(nonactset_len, 0, sqrt(sigmaalpha));
//   if(nonactset_len > 0 ){
//     norm = rnorm(nonactset_len, 0, sqrt(sigmaalpha));
//     alpha.elem(nonactset) = norm;
//   }
//   if(actset_len>0){
//     if(actset_len > 1) {
//       mat XA = X.cols(actset);
//       mat premat = (XA.t()*XA)/sigmae;
//       premat.diag() = premat.diag() + (1/sigmaalpha);
//       mat premat_chol = chol(premat);
//       vec mu = (XA.t()*y)/sigmae;
//       vec b = solve(trimatu(premat_chol).t(), mu);
//       //Rcout << "b" << b;
//       vec Z = rnorm(actset_len, 0, 1);
//      // Rcout << "Z" << Z;
//       alpha.elem(actset) = solve(trimatu(premat_chol), Z+b);
//       return(alpha);
//     } else {
//       vec XA = X.col(actset(0));
//       vec XAsq = pow(XA,2);
//       double XAsum = sum(XAsq);
//       double varx = 1/(XAsum + 1/sigmaalpha);
//       double meanx = varx*sum(XA.t()*y)/sigmae;
//       norm = rnorm(1,meanx, sqrt(varx));
//       alpha.elem(actset) = norm;
//       return(alpha);
//     }
//   }
//   return(alpha);
//   
// }
// 
// //arma::vec update_sigma(vec curr.gamma, vec lambda, arma::mat X) {
//   //return(sigma);
// //}
// 
// // [[Rcpp::export]]
// List update_epsilon(double epsilon,  double tauepsilon, double sigmagamma, double esd, double emu,
//                          arma::vec gamma, arma::vec eigenvalue,
//                          arma::mat laplacian, arma::mat Ip) {
//   double epsilonnew = rnorm(1, epsilon, sqrt(tauepsilon))(0);
//   vec sqgamma = pow(gamma,2);
//   double sgamma = sum(sqgamma);
//   double elikecurr = 0.5*sum(log(eigenvalue+exp(epsilon))) - exp(epsilon)*sgamma*0.5/sigmagamma - epsilon - 0.5*pow(epsilon-emu, 2)/(esd*esd);
//   mat eigmatnew = laplacian + exp(epsilonnew)*Ip;
//   double elikenew = 0.5*sum(log(eigenvalue+exp(epsilonnew))) - exp(epsilonnew)*sgamma*0.5/sigmagamma-epsilonnew- 0.5*pow(epsilonnew-emu, 2)/(esd*esd) ;
//   double ediff = elikenew - elikecurr;
//   List L =List::create(Named("ediff")=ediff,
//                        Named("epsilonnew")= epsilonnew,
//                        Named("eigmatnew")= eigmatnew);
//   return(L);
// }
// 
// //arma::vec update_lambda(vec currgamma, vec lambda, arma::mat X) {
//   //dtruncnorm(1, 2,3,4,5);
//   //return(sigma);
// //}