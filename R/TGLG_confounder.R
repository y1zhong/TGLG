TGLG_confounder = function(X, cf, y, y.type='continuous', net=NULL,nsim=30000, ntune=10000, freqTune=100,
                                   nest=10000, nthin=10, a.alpha=0.01,b.alpha=0.01,
                                   a.gamma=0.01,b.gamma=0.01,
                                   ae0=0.01,be0=0.01,lambda.lower=0,
                                   lambda.upper=10, emu=-5, esd=3,prob_select=0.95,seed=123){
  if(nrow(X) != nrow(cf)) stop('rownames of X and cf do not match')
  X_full <- as.matrix(cbind(X,cf))
  net_full <- net + vertices(colnames(cf))
  if(y.type == 'binary') {
    return(TGLG_binary(X_full, y, net_full, nsim, ntune, freqTune,nest,
                       nthin, a.alpha,b.alpha,a.gamma=,b.gamma,
                       lambda.lower, lambda.upper, emu, esd,prob_select,seed))
  } else {
    return(TGLG_continuous_revised_sampling(X=X_full, y=y, net=net_full, nsim, ntune, freqTune,nest,
                       nthin, a.alpha,b.alpha,a.gamma,b.gamma,
                       lambda.lower, lambda.upper, emu, esd,prob_select,seed))
  }
}