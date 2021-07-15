samp_mvn <- function(mu, sigma){
  U <- chol(sigma)
  X <- mu+crossprod(U)
}