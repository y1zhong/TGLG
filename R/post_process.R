#' @export
prob_fdr = function(tglg, rate = 0.1){
  post <- tglg$post_summary
  prob <- post$selectProb
  sort.p <- sort(prob,decreasing = T)
  lambda.vec <- sapply(1:length(prob), function(x) sum(1-sort.p[1:x])/x )
  lambda <- max(which(lambda.vec<rate))
  phi_alpha <- sort.p[lambda]
  return(phi_alpha)
}


inclusion_prob = function(tglg, v){
  
}