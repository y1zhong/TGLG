# TGLG model fitting for continous outcome
# Input Arguments:
#   X: input features, dimension n*p, with n as sample size and p as number of features
#   y: response variable: continous value
#   net: igraph object that represents the network
#   nsim: total number of MCMC iteration
#   ntune: first number of MCMC iteration to adjust propose variance to get desired acceptance rate
#   freqTune: frequency to tune the variance term.
#   nest: number used for estimation
#   nthin: thin MCMC chain
#   a.alpha, b.alpha: inverse-gamma distribution parameter for sigma_alpha
#   a.gamma, b.gamma: inverse-gamma distribution parameter for sigma_gamma
#   lambda.lower, lambda.uppper: lower and upper bound of uniform distribution for lambda
#   emu, esd: mean and standard deviation of log-normal distribution for epsilon
#
# Output object:
#   post_summary: a dataframe including the selection probablity and posterior mean of betas
#   dat: X, y, and net used to generate results
#   save_mcmc: all the mcmc samples saved after burnin and thinning.
#' @export
#' @import glmnet
#' @import mvnfast
#' @import MCMCpack
#' @import pROC
#' @import truncnorm
#' @import Matrix
TGLG_continuous_conf = function(X, y, net=NULL,nsim=30000, ntune=10000, freqTune=100,
                           nest=10000, nthin=10,  a.alpha=0.01,b.alpha=0.01,
                           a.gamma=0.01,b.gamma=0.01,
                           ae0=0.01,be0=0.01,lambda.lower=0,
                           lambda.upper=10, emu=-5, esd=3,prob_select=0.95,
                           seed=123, confound, beta.ini.meth='LASSO',
                           lambda.ini='median', noise=FALSE, select.prop=0.9){
  # if(X.scale == T) X <- scale(X)
  # if(!is.null(confound)){
  #   if (sum(apply(confound, 2, class) != "numeric") != 0) {
  #     stop('confounder class type not numeric')
  #   }
  #   if(nrow(X) != nrow(confound)) stop('rownames of X and cf do not match')
  #   X_feat <- X
  #   #confound_scale <- scale(confound)
  #   p_feat <- ncol(X_feat)
  #   X <- as.matrix(cbind(X,confound))
  #   net <- net + vertices(colnames(confound))
  # }
  set.seed(seed)
  p = ncol(X) #number of features
  if(is.null(net)){
    adjmat = diag(p)
    net=graph.adjacency(adjmat,mode="undirected")
    
  }
  laplacian = laplacian_matrix(net,normalized=TRUE,sparse=FALSE)
  #laplacian_sp <- Matrix(laplacian, sparse = TRUE)
  
  tau.gamma=0.01/p
  tau.lambda=0.1
  tau.epsilon = 0.5
  tau.omega = 20
  #count number of acceptance during MCMC updates
  accept.gamma=0
  accept.epsilon=0
  accept.lambda =0
  accept.omega =0
  #get initial value using lasso
  #alpha = 0.5
  if(beta.ini.meth == 'LASSO'){
    beta.elas=cv.glmnet(X,y,alpha=1,intercept=FALSE)
    beta.ini=as.numeric(coef(beta.elas,beta.elas$lambda.min))[-1]
  } else if(beta.ini.meth  == 'Univariate') {
    #inital beta values for features
    beta.ini.feat <- sapply(1:p, function(k){
      beta_var <- cbind(X[, k], confound)
      fit_temp <- glm(y~., data=data.frame(beta_var))
      unname(fit_temp$coefficients[2])
    })

    #fit_temp <- glm(y~., data=data.frame(confound))
    #beta.ini.conf <-  unname(fit_temp$coefficients[-1])
    # beta.ini <- c(beta.ini.feat, beta.ini.conf)
    beta.ini <- beta.ini.feat
  } else if(beta.ini.meth  == 'Mixed') {
    #inital beta values for features
    beta.ini.feat <- sapply(1:p_feat, function(k){
      beta_var <- cbind(X_feat[, k], confound)
      fit_temp <- glm(y~., data=data.frame(beta_var))
      unname(fit_temp$coefficients[2])
    })
    
  } else {
    stop('beta ini method input incorrect')
  }
  
  if(noise == T){
    set.seed(seed)
    noise.add <- rnorm(p, 0, 0.5)
    beta.ini = beta.ini+noise.add
  }
  
  
  #initial value
  sigmae=sd(y)#rinvgamma(1,ae0,be0)
  betaini2 = beta.ini[which(beta.ini!=0)]
  if(length(betaini2)==0){
    betaini2=0.001
  }
  #change lambda inital
  if(lambda.ini == 'median'){
    lambda =  median(abs(betaini2))/2
  } else if (lambda.ini == 'zero') {
    lambda = 0
  } else if (lambda.ini == 'prop') {
    beta.abs <- abs(beta.ini)
    lambda = quantile(beta.abs, select.prop)
  } else {
    stop('lambda initial incorrect')
  }
  
  burnin <- nsim-nest
  gamma =beta.ini
  alpha =beta.ini
  beta = alpha*as.numeric(abs(gamma)>lambda)
  sigma.alpha = 5
  sigma.gamma = b.gamma/a.gamma
  epsilon = emu
  eigen.value = eigen_val(laplacian)
  Ip = diag(1,nrow(laplacian)) #identy matrix

  eigmat = laplacian+exp(epsilon)*Ip #inverse of graph laplacian matrix

  #confounders
  sigma.omega <- 50
  Z <- as.matrix(confound)
  p.conf <- ncol(Z)
  prematZ = crossprod(Z)/sigmae
  cholZ <- chol(prematZ)
  diag(prematZ) <- diag(prematZ) + (1/sigma.omega)
  diag.prematZ <-diag(prematZ)
  omega <- rnorm(p.conf, 0, sd=sqrt(sigma.omega))
  
  #matrix to store values for all iterations
  #nmcmc <- nest/nthin
  numrow = nest/nthin
  fgamma=matrix(0,nrow=numrow,ncol=length(gamma))
  
  Beta =matrix(0,nrow=numrow,ncol=p)
  Alpha =matrix(0,nrow=numrow,ncol=p)
  Gamma = matrix(0,nrow=numrow,ncol=p)
  Lambda=rep(0,numrow)
  Sigma.gamma=rep(0,numrow)
  Sigma.alpha = rep(0, numrow)
  Epsilon = rep(0, numrow)
  likelihood = rep(0, numrow)
  Sigmae = rep(0, numrow)
  #for confounders
  Omega = matrix(0,nrow=numrow,ncol=p.conf)
  
  iter = seq((nsim-nest+1),nsim,by=nthin)
  #SimTime = rep(0, nsim)
  #sim=1
  pb = txtProgressBar(style=3)
  for(sim in 1:nsim){
    #y_temp <- y
    
    #update gamma using MALA
    curr.gamma=gamma
    curr.hgamma= curr.gamma^2-lambda^2
    curr.beta=alpha*as.numeric(abs(curr.gamma)>lambda)
    
    
    #update omega
    yz <- y - X%*%curr.beta
    b <- backsolve(premat_chol, mu, transpose=T)
    r <- rnorm(actset_len, 0, 1)
    alpha[actset] <- backsolve(premat_chol, r+b)
   # for(i in 1:p.conf){
   #    z <- rnorm(1, 0, 1)
   #    aii <- diag.prematZ[i]
   #    if(i == 1) omega[i] <- z/sqrt(aii) -sum((prematZ[i, (i+1):p.conf]*omega[(i+1):p.conf]))/aii; next
   #    if(i == p.conf) omega[i] <- z/sqrt(aii) -sum((-prematZ[i, 1:(i-1)]*omega[1:(i-1)]))/aii; next
   #    omega[i] <- z/sqrt(aii) -(sum((prematZ[i, (i+1):p.conf]*omega[(i+1):p.conf])) - 
   #      sum(prematZ[i, 1:(i-1)]*omega[1:(i-1)]))/aii
   #  }
    
    # y <- y - X%*%curr.beta
    # yerr.omega <- y - Z %*% omega
    # new.omega <- curr.omega+rnorm(length(omega),0, sd=sqrt(tau.omega))
    # #curr.adiff.omega = new.omega - curr.omega
    # curr.den.omega = 0.5*sum(curr.omega^2)/sigma.omega
    # curr.post.omega=0.5*sum(yerr.omega^2)/sigmae
    # 
    # new.yerr.omega <- y - Z %*% new.omega
    # new.den.omega = 0.5*sum(new.omega^2)/sigma.omega
    # new.post.omega=0.5*sum(new.yerr.omega^2)/sigmae
    # u=runif(1)
    # post.diff.omega=new.post.omega+new.den.omega-curr.post.omega-curr.den.omega
    # #cat(paste('omega',post.diff.omega),'\n')
    # if(post.diff.omega>=log(u)){
    #   accept.omega=accept.omega+1
    #   omega=new.omega
    # }
    
    y <- y_temp
    y <- y-Z%*%omega
    yerr = y-X%*%curr.beta
    curr.diff=crossprod(X, yerr)
    curr.dgamma = dappro2(curr.hgamma)*curr.gamma
    
    tau.gamma_half <- tau.gamma/2
    sigmae_twice <- sigmae*2
    
    eigmat_curr.gamma <- eigmat%*%curr.gamma
    #curr.post=crossprod(yerr)/sigmae*0.5+0.5*crossprod(curr.gamma, crossprod(eigmat, curr.gamma))/sigma.gamma
    curr.post=sum(yerr^2)/sigmae_twice+0.5*crossprod(curr.gamma, eigmat_curr.gamma)/sigma.gamma
    new.gamma.mean = curr.gamma + tau.gamma_half*(curr.diff*gamma*curr.dgamma/sigmae - eigmat_curr.gamma/sigma.gamma)
    new.gamma = new.gamma.mean+rnorm(length(gamma),0, sd=sqrt(tau.gamma))
    new.gamma = as.numeric(new.gamma)
    new.hgamma = new.gamma^2-lambda^2
    new.beta=alpha*as.numeric(abs(new.gamma)>lambda)
    new.yerr = y-X%*%new.beta
    new.diff = crossprod(X, new.yerr)
    new.dgamma = dappro2(new.hgamma)*new.gamma
    curr.gamma.mean = new.gamma + tau.gamma_half*(new.diff*gamma*new.dgamma/sigmae - eigmat%*%new.gamma/sigma.gamma)
    
    #new.post=crossprod(new.yerr)/sigmae*0.5+0.5*crossprod(new.gamma, crossprod(eigmat, new.gamma))/sigma.gamma
    new.post=sum(new.yerr^2)/sigmae_twice+0.5*crossprod(new.gamma, crossprod(eigmat, new.gamma))/sigma.gamma
    curr.adiff = new.gamma- new.gamma.mean
    curr.den = 0.5*sum(curr.adiff^2)/tau.gamma
    new.adiff=curr.gamma-curr.gamma.mean
    new.den =0.5* sum(new.adiff^2)/tau.gamma
    u=runif(1)
    post.diff=curr.post+curr.den-new.post-new.den
    #cat(paste('gamma',post.diff[1,1]),'\n')
    if(post.diff>=log(u)){
      accept.gamma=accept.gamma+1
      gamma=new.gamma
    }
    beta=alpha*as.numeric(abs(gamma)>lambda)
    #update alpha
    
    actset=which(abs(gamma)>lambda)
    non.actset=setdiff(1:p,actset)
    actset_len <- length(actset)
    if(length(non.actset)>0 ){
      alpha[non.actset]=rnorm(length(non.actset),mean=0,sd=sqrt(sigma.alpha))
    }
    
    if(actset_len>0){
      if(actset_len > 1) {
        XA=X[,actset]
        premat = crossprod(XA)/sigmae
        diag(premat) <- diag(premat) + (1/sigma.alpha)
        
        premat_chol <- chol(premat)
        mu <- crossprod(XA, y)/sigmae
        #transpose in function
        #b2 <- forwardsolve(premat_chol, mu, transpose = T)
        b <- backsolve(premat_chol, mu, transpose=T)
        r <- rnorm(actset_len, 0, 1)
        alpha[actset] <- backsolve(premat_chol, r+b)
        
      } else {
        XA = X[, actset]
        varx = 1/(sum(XA^2) + 1/sigma.alpha)
        meanx = varx*sum(XA*y)/sigmae
        alpha[actset] = rnorm(1, mean = meanx, sd = sqrt(varx))
      }
    }
    
    #update sigma.gamma
    a.gamma.posterior=a.gamma+0.5*p
    b.gamma.posterior=b.gamma+0.5*crossprod(gamma, eigmat%*%gamma)
    sigma.gamma=rinvgamma(1,a.gamma.posterior,b.gamma.posterior)
    
    beta=alpha*as.numeric(abs(gamma)>lambda)
    a.alpha.posterior = a.alpha + 0.5*p
    b.alpha.posterior = b.alpha + 0.5*sum(alpha^2)
    sigma.alpha = rinvgamma(1, a.alpha.posterior, a.alpha.posterior)
    
    ae=ae0+nrow(X)/2
    be=be0+0.5*sum((y-X%*%beta)^2)
    sigmae=rinvgamma(1,ae,be)
    
    #update epsilon
    epsilon.new = rnorm(1,mean = epsilon, sd = sqrt(tau.epsilon))
    sgamma = sum(gamma^2)
    
    sgamma_cal <- sgamma*0.5/sigma.gamma
    esd_square <- esd^2
    
    elike.curr = 0.5*sum(log(eigen.value+exp(epsilon))) - exp(epsilon)*sgamma_cal - epsilon - 0.5*(epsilon-emu)^2/esd_square
    eigmat.new = laplacian + exp(epsilon.new)*Ip
    elike.new = 0.5*sum(log(eigen.value+exp(epsilon.new))) - exp(epsilon.new)*sgamma_cal-epsilon.new- 0.5*(epsilon.new-emu)^2/esd_square
    eu = runif(1)
    ediff = elike.new - elike.curr
    if(ediff>=log(eu)) {
      epsilon = epsilon.new
      eigmat = eigmat.new
      accept.epsilon = accept.epsilon+1
    }
    
    #update lambda
    lambda.new=rtruncnorm(1,a=lambda.lower,b=lambda.upper,mean=lambda,sd=sqrt(tau.lambda))
    #curr.lpost=-crossprod(y-X%*%beta)/sigmae*0.5
    curr.lpost=-sum((y-X%*%beta)^2)/sigmae*0.5
    new.beta=alpha*as.numeric(abs(gamma)>lambda.new)
    new.lpost=-sum((y-X%*%new.beta)^2)/sigmae*0.5
    dnew = log(dtruncnorm(lambda.new,a=lambda.lower,b=lambda.upper,mean=lambda,sd=sqrt(tau.lambda)))
    dold = log(dtruncnorm(lambda,a=lambda.lower,b=lambda.upper,mean=lambda.new,sd=sqrt(tau.lambda)))
    llu=runif(1)
    ldiff=new.lpost + dold-curr.lpost-dnew
    if(ldiff>log(llu)){
      lambda=lambda.new
      beta=new.beta
      accept.lambda = accept.lambda+1
    }
    
    if(sim<=ntune&sim%%freqTune==0){
      tau.gamma=adjust_acceptance(accept.gamma/100,tau.gamma,0.5)
      tau.lambda=adjust_acceptance(accept.lambda/100,tau.lambda,0.3)
      tau.epsilon=adjust_acceptance(accept.epsilon/100,tau.epsilon,0.3)
      tau.omega=adjust_acceptance(accept.omega/100,tau.omega,0.5)
      accept.gamma=0
      accept.epsilon=0
      accept.lambda =0
      
    }
    if(sim %in% iter){
      idx <- match(sim, iter)
      Gamma[idx,]=gamma
      Alpha[idx, ] = alpha
      Sigma.gamma[idx] = sigma.gamma
      Sigma.alpha[idx] = sigma.alpha
      Sigmae[idx] =sigmae
      Epsilon[idx] = epsilon
      Lambda[idx]=lambda
      Beta[idx,]=beta
      Omega[idx,]=omega
      likelihood[idx] = -sum((y-X%*%beta)^2)/sigmae*0.5
    }
    #change from 0 to 1 to match iter index
    if(sim>burnin & ((sim-burnin) %% nthin == 1)){
      j=sim-burnin
      fgamma[j/nthin,]=as.numeric(abs(gamma)>lambda)
    }
    y <- y_temp
   # print(y[1])
    setTxtProgressBar(pb,sim/nsim)
  }
  
  close(pb)
  #use last nest iterations to do the estimation and thin by nthin
  
  #selection probability for each predictor
  gammaSelProb=apply(fgamma,2,mean)
  
  hbeta=Beta
  
  gammaSelId=as.numeric(gammaSelProb>prob_select)
  beta.est=rep(0,p)
  beta.select=which(gammaSelId==1)
  for(j in 1:length(beta.select)) {
    colid =beta.select[j]
    beta.est[colid]=mean(hbeta[fgamma[,colid]==1, colid])
  }
  
  post_summary = data.frame(selectProb = gammaSelProb, betaEst = beta.est)
  #iter = seq((nsim-nest+1),nsim,by=nthin)
  save_mcmc = cbind(iter,Beta,Alpha,Gamma,Omega,
                    Lambda,Sigma.gamma,Sigma.alpha,Epsilon,Sigmae,
                    likelihood)
  colnames(save_mcmc) = c("iter",
                          paste("beta",1:p,sep=""),
                          paste("alpha",1:p,sep=""),
                          paste("gamma",1:p,sep=""),
                          paste("omega",1:p.conf,sep=""),
                          "lambda","sigma_gamma","sigma_alpha","epsilon","sigmae","loglik")
  
  return(list(post_summary=post_summary, dat = list(X=X,y=y,net=net), save_mcmc = save_mcmc))
}
