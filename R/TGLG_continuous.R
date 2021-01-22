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
TGLG_continuous = function(X, y, net=NULL,nsim=30000, ntune=10000, freqTune=100,
                           nest=10000, nthin=10, a.alpha=0.01,b.alpha=0.01,
                           a.gamma=0.01,b.gamma=0.01,
                           ae0=0.01,be0=0.01,lambda.lower=0,
                           lambda.upper=10, emu=-5, esd=3,prob_select=0.95){

  p = ncol(X) #number of features

  if(is.null(net)){
    adjmat = diag(p)
    net=graph.adjacency(adjmat,mode="undirected")

  }
  laplacian = laplacian_matrix(net,normalized=TRUE,sparse=FALSE)

  tau.gamma=0.01/p
  tau.lambda=0.1
  tau.epsilon = 0.5
  #count number of acceptance during MCMC updates
  accept.gamma=0
  accept.epsilon=0
  accept.lambda =0
  #get initial value using lasso
  beta.elas=cv.glmnet(X,y,alpha=1,intercept=FALSE)
  beta.ini=as.numeric(coef(beta.elas,beta.elas$lambda.min))[-1]
  #initial value
  sigmae=sd(y)#rinvgamma(1,ae0,be0)
  betaini2 = beta.ini[which(beta.ini!=0)]
  if(length(betaini2)==0){
    betaini2=0.001
  }
  lambda =  median(abs(betaini2))/2
  gamma =beta.ini
  alpha =beta.ini
  beta = alpha*as.numeric(abs(gamma)>lambda)
  sigma.alpha = 5
  sigma.gamma = b.gamma/a.gamma
  epsilon = emu
  eigen.value = eigen(laplacian)$values
  Ip = diag(1,nrow(laplacian)) #identy matrix

  eigmat = laplacian+exp(epsilon)*Ip #inverse of graph laplacian matrix

  #matrix to store values for all iterations
  Beta =matrix(0,nrow=nsim,ncol=p)
  Alpha =matrix(0,nrow=nsim,ncol=p)
  Gamma = matrix(0,nrow=nsim,ncol=p)
  Lambda=rep(0,nsim)
  Sigma.gamma=rep(0,nsim)
  Sigma.alpha = rep(0, nsim)
  Epsilon = rep(0, nsim)
  likelihood = rep(0, nsim)
  Sigmae = rep(0, nsim)
  pb = txtProgressBar(style=3)
  for(sim in 1:nsim){
    #update gamma using MALA
    curr.gamma=gamma
    curr.hgamma= curr.gamma^2-lambda^2
    curr.beta=alpha*as.numeric(abs(curr.gamma)>lambda)
    yerr = y-X%*%curr.beta
    curr.diff=crossprod(X, yerr)
    curr.dgamma = dappro2(curr.hgamma)*curr.gamma
    curr.post=crossprod(yerr)/sigmae*0.5+0.5*crossprod(curr.gamma, crossprod(eigmat, curr.gamma))/sigma.gamma
    new.gamma.mean =curr.gamma + tau.gamma/2*(curr.diff*gamma*curr.dgamma/sigmae - eigmat%*%curr.gamma/sigma.gamma)
    new.gamma = new.gamma.mean+rnorm(length(gamma),0, sd=sqrt(tau.gamma))
    new.gamma = as.numeric(new.gamma)
    new.hgamma = new.gamma^2-lambda^2
    new.beta=alpha*as.numeric(abs(new.gamma)>lambda)
    new.yerr = y-X%*%new.beta
    new.diff = crossprod(X, new.yerr)
    new.dgamma = dappro2(new.hgamma)*new.gamma
    curr.gamma.mean = new.gamma + tau.gamma/2*(new.diff*gamma*new.dgamma/sigmae - eigmat%*%new.gamma/sigma.gamma)
    new.post=crossprod(new.yerr)/sigmae*0.5+0.5*crossprod(new.gamma, crossprod(eigmat, new.gamma))/sigma.gamma
    curr.adiff = new.gamma- new.gamma.mean
    curr.den = 0.5*sum(curr.adiff^2)/tau.gamma
    new.adiff=curr.gamma-curr.gamma.mean
    new.den =0.5* sum(new.adiff^2)/tau.gamma
    u=runif(1)
    post.diff=curr.post+curr.den-new.post-new.den
    if(post.diff>=log(u)){
      accept.gamma=accept.gamma+1
      gamma=new.gamma
    }
    Gamma[sim,]=gamma
    beta=alpha*as.numeric(abs(gamma)>lambda)
    #update alpha

    actset=which(abs(gamma)>lambda)
    non.actset=setdiff(1:p,actset)

    if(length(non.actset)>0 ){
      alpha[non.actset]=rnorm(length(non.actset),mean=0,sd=sqrt(sigma.alpha))
    }

    if(length(actset)>0){
      if(length(actset) > 1) {
        XA=X[,actset]
        pvar = diag(rep(1/sigma.alpha,length(actset)))
        premat = crossprod(XA)/sigmae + pvar
        varx = chol2inv(chol(premat))
        meanx = varx%*%crossprod(XA, y)/sigmae
        alpha[actset] = rmvn(1,  mu=meanx, sigma = varx)
        #gamma[actset] = rmvnorm(1, mean=meanx, sigma = varx)
      } else {
        XA = X[, actset]
        varx = 1/(sum(XA^2) + 1/sigma.alpha)
        meanx = varx*sum(XA*y)/sigmae
        alpha[actset] = rnorm(1, mean = meanx, sd = sqrt(varx))
      }
    }
    Alpha[sim, ] = alpha
    #update sigma.gamma
    a.gamma.posterior=a.gamma+0.5*p
    b.gamma.posterior=b.gamma+0.5*crossprod(gamma,crossprod(eigmat,gamma))
    sigma.gamma=rinvgamma(1,a.gamma.posterior,b.gamma.posterior)
    Sigma.gamma[sim] = sigma.gamma

    beta=alpha*as.numeric(abs(gamma)>lambda)
    a.alpha.posterior = a.alpha + 0.5*p
    b.alpha.posterior = b.alpha + 0.5*sum(alpha^2)
    sigma.alpha = rinvgamma(1, a.alpha.posterior, a.alpha.posterior)
    Sigma.alpha[sim] = sigma.alpha

    ae=ae0+nrow(X)/2
    be=be0+0.5*crossprod(y-X%*%beta)
    sigmae=rinvgamma(1,ae,be)
    Sigmae[sim] =sigmae

    #update epsilon
    epsilon.new = rnorm(1,mean = epsilon, sd = sqrt(tau.epsilon))
    sgamma = sum(gamma^2)
    elike.curr = 0.5*sum(log(eigen.value+exp(epsilon))) - exp(epsilon)*sgamma*0.5/sigma.gamma - epsilon - 0.5*(epsilon-emu)^2/esd^2
    eigmat.new = laplacian + exp(epsilon.new)*Ip
    elike.new = 0.5*sum(log(eigen.value+exp(epsilon.new))) - exp(epsilon.new)*sgamma*0.5/sigma.gamma-epsilon.new- 0.5*(epsilon.new-emu)^2/esd^2
    eu = runif(1)
    ediff = elike.new - elike.curr
    if(ediff>=log(eu)) {
      epsilon = epsilon.new
      eigmat = eigmat.new
      accept.epsilon = accept.epsilon+1
    }
    Epsilon[sim] = epsilon

    #update lambda
    lambda.new=rtruncnorm(1,a=lambda.lower,b=lambda.upper,mean=lambda,sd=sqrt(tau.lambda))
    curr.lpost=-crossprod(y-X%*%beta)/sigmae*0.5
    new.beta=alpha*as.numeric(abs(gamma)>lambda.new)
    new.lpost=-crossprod(y-X%*%new.beta)/sigmae*0.5
    dnew = log(dtruncnorm(lambda.new,a=lambda.lower,b=lambda.upper,mean=lambda,sd=sqrt(tau.lambda)))
    dold = log(dtruncnorm(lambda,a=lambda.lower,b=lambda.upper,mean=lambda.new,sd=sqrt(tau.lambda)))
    lu=runif(1)
    ldiff=new.lpost + dold-curr.lpost-dnew
    if(ldiff>log(lu)){
      lambda=lambda.new
      beta=new.beta
      accept.lambda = accept.lambda+1
    }
    Lambda[sim]=lambda
    Beta[sim,]=beta

    likelihood[sim] = -crossprod(y-X%*%beta)/sigmae*0.5


    if(sim<=ntune&sim%%freqTune==0){
      tau.gamma=adjust_acceptance(accept.gamma/100,tau.gamma,0.5)
      tau.lambda=adjust_acceptance(accept.lambda/100,tau.lambda,0.3)
      tau.epsilon=adjust_acceptance(accept.epsilon/100,tau.epsilon,0.3)
      accept.gamma=0
      accept.epsilon=0
      accept.lambda =0

    }

    setTxtProgressBar(pb,sim/nsim)
  }

  close(pb)
  #use last nest iterations to do the estimation and thin by nthin
  numrow = nest/nthin
  fgamma=matrix(0,nrow=numrow,ncol=length(gamma))
  for(j in 1:nest){
    if(j%%nthin==0){
      fgamma[j/nthin,]=as.numeric(abs(Gamma[nsim-nest+j,])>Lambda[nsim-nest+j])
    }
  }
  #selection probability for each predictor
  gammaSelProb=apply(fgamma,2,mean)

  hbeta=Beta[((nsim-nest+1):nsim)%%nthin==1,]

  gammaSelId=as.numeric(gammaSelProb>prob_select)
  beta.est=rep(0,p)
  beta.select=which(gammaSelId==1)
  for(j in 1:length(beta.select)) {
    colid =beta.select[j]
    beta.est[colid]=mean(hbeta[fgamma[,colid]==1, colid])
  }

  post_summary = data.frame(selectProb = gammaSelProb, betaEst = beta.est)
  iter = seq((nsim-nest+1),nsim,by=nthin)
  save_mcmc = cbind(iter,Beta[iter,],Alpha[iter,],Gamma[iter,],
                    Lambda[iter],Sigma.gamma[iter],Sigma.alpha[iter],Epsilon[iter],Sigmae[iter],
                    likelihood[iter])
  colnames(save_mcmc) = c("iter",
                          paste("beta",1:p,sep=""),
                          paste("alpha",1:p,sep=""),
                          paste("gamma",1:p,sep=""),
                          "lambda","sigma_gamma","sigma_alpha","epsilon","sigmae","loglik")

  return(list(post_summary=post_summary, dat = list(X=X,y=y,net=net), save_mcmc = save_mcmc))
}
