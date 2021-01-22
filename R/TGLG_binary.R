#' TGLG model fitting for binary outcome
#'
#' @param X: input features, dimension n*p, with n as sample size and p as number of features
#' @param y: response variable: binary 0, 1
#' @param net: igraph object that represents the network
#' @param nsim: total number of MCMC iteration
#' @param ntune: first number of MCMC iteration to adjust propose variance to get desired acceptance rate
#' @param freqTune: frequency to tune the variance term.
#' @param nest: number used for estimation
#' @param nthin: thin MCMC chain
#' @param a.alpha, b.alpha: inverse-gamma distribution parameter for sigma_alpha
#' @param a.gamma, b.gamma: inverse-gamma distribution parameter for sigma_gamma
#' @param lambda.lower, lambda.uppper: lower and upper bound of uniform distribution for lambda
#' @param emu, esd: mean and standard deviation of log-normal distribution for epsilon
#' @return post_summary: a dataframe including the selection probablity and posterior mean of betas
#' @return dat: X, y, and net used to generate results
#' @return save_mcmc: all the mcmc samples saved after burnin and thinning.
#' @export
#'
TGLG_binary = function(X, y, net=NULL,nsim=30000, ntune=10000, freqTune=100,nest=10000,
                       nthin=10, a.alpha=0.01,b.alpha=0.01,a.gamma=0.01,b.gamma=0.01,
                       lambda.lower=0, lambda.upper=10, emu=-5, esd=3,prob_select=0.95){
  p = ncol(X) #number of features
  if(is.null(net)){
    adjmat = matrix(0,nrow=p,ncol=p)
    net=graph.adjacency(adjmat,mode="undirected")

  }

  laplacian = laplacian_matrix(net,normalized=TRUE,sparse=FALSE)
  #proposed variance for
  tau.gamma=0.01/p
  tau.alpha = 0.1/sqrt(p)
  tau.lambda=0.1
  tau.epsilon = 0.5
  #count number of acceptance during MCMC updates
  accept.gamma=0
  accept.alpha=0
  accept.epsilon=0
  accept.lambda =0
  #get initial value using lasso
  #type.measure = "deviance"
  #if(family=="binomial"){
  #type.measure="class"
  #}
  beta.elas=cv.glmnet(X,y,family="binomial",alpha=1,intercept=FALSE,type.measure="class")
  beta.ini=as.numeric(coef(beta.elas,beta.elas$lambda.min))[-1]
  #initial value
  betaini2 = beta.ini[which(beta.ini!=0)]
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
  pb = txtProgressBar(style = 3)
  for(sim in 1:nsim){
    #update gamma using MALA
    curr.gamma=gamma
    curr.hgamma= curr.gamma^2-lambda^2
    curr.beta=alpha*as.numeric(abs(curr.gamma)>lambda)
    curr.dgamma = dappro2(curr.hgamma)*curr.gamma
    curr.temp = X%*%curr.beta
    curr.post=sum(y*curr.temp) - sum(log(1+exp(curr.temp))) - 0.5*curr.gamma%*%(eigmat%*%curr.gamma)/sigma.gamma
    currz = 1/(1+exp(-curr.temp))
    curr.diff = crossprod(X, y-currz)
    new.gamma.mean = curr.gamma + tau.gamma/2*(curr.diff*gamma*curr.dgamma - crossprod(eigmat,curr.beta)/sigma.gamma)
    new.gamma.mean =as.numeric(new.gamma.mean)
    new.gamma=new.gamma.mean + rnorm(p, mean=0, sd =sqrt(tau.gamma))
    new.beta=alpha*as.numeric(abs(new.gamma)>lambda)
    new.temp = X%*%new.beta
    new.hgamma= new.gamma^2-lambda^2
    #new.beta=gamma*as.numeric(abs(new.gamma)>lambda)
    new.dgamma = dappro2(new.hgamma)*new.gamma
    newz = 1/(1+exp(-new.temp))
    new.diff = crossprod(X, y-newz)
    curr.gamma.mean = new.gamma + tau.gamma/2*(new.diff*gamma*new.dgamma - crossprod(eigmat,new.beta)/sigma.gamma)
    new.post=sum(y*new.temp) - sum(log(1+exp(new.temp))) - 0.5*new.gamma%*%(eigmat%*%new.gamma)/sigma.gamma
    u=runif(1)
    post.diff=new.post - 0.5*sum((curr.gamma-curr.gamma.mean)^2/tau.gamma)-curr.post +
      0.5*sum((new.gamma-new.gamma.mean)^2/tau.gamma)
    if(post.diff>=log(u)){
      accept.gamma=accept.gamma+1
      gamma = new.gamma
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
      XA = X[, actset,drop=F]
      calpha = alpha[actset]
      gcurr.temp = XA%*%calpha
      curr.gpost=sum(y*gcurr.temp) - sum(log(1+exp(gcurr.temp))) - 0.5*sum(calpha^2)/sigma.alpha
      gcurrz = 1/(1+exp(-gcurr.temp))
      new.alpha.mean =calpha + tau.alpha/2*(crossprod(XA, y-gcurrz) -calpha/sigma.alpha )
      new.alpha.mean = as.numeric(new.alpha.mean)
      new.aalpha = new.alpha.mean + rnorm(length(actset), mean =0, sd =sqrt(tau.alpha))
      gnew.temp = XA%*%new.aalpha
      new.gpost=sum(y*gnew.temp) - sum(log(1+exp(gnew.temp))) - 0.5*sum(new.aalpha^2)/sigma.alpha
      gnewz = 1/(1+exp(-gnew.temp))
      curr.alpha.mean =new.aalpha + tau.alpha/2*(crossprod(XA, y-gnewz) -new.aalpha/sigma.alpha )
      gu=runif(1)
      gpost.diff=new.gpost- 0.5*sum((calpha-curr.alpha.mean)^2/tau.alpha)-curr.gpost +
        0.5*sum((new.aalpha-new.alpha.mean)^2/tau.alpha)
      if(gpost.diff>=log(gu)){
        accept.alpha=accept.alpha+1
        alpha[actset] = new.aalpha
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

    lambda.new=rtruncnorm(1,a=lambda.lower,b=lambda.upper,mean=lambda,sd=sqrt(tau.lambda))
    curr.ltemp = X%*%beta
    curr.lpost=sum(y*curr.ltemp) - sum(log(1+exp(curr.ltemp)))
    new.beta=alpha*as.numeric(abs(gamma)>lambda.new)
    new.ltemp = X%*%new.beta
    new.lpost=sum(y*new.ltemp) - sum(log(1+exp(new.ltemp)))
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

    temp = X%*%beta
    likelihood[sim]= sum(y*temp) - sum(log(1+exp(temp)))


    if(sim<=ntune&sim%%freqTune==0){
      tau.gamma=adjust_acceptance(accept.gamma/100,tau.gamma,0.5)
      tau.lambda=adjust_acceptance(accept.lambda/100,tau.lambda,0.3)
      tau.alpha=adjust_acceptance(accept.alpha/100,tau.alpha,0.5)
      tau.epsilon=adjust_acceptance(accept.epsilon/100,tau.epsilon,0.3)
      accept.gamma=0
      accept.alpha=0
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
                    Lambda[iter],Sigma.gamma[iter],Sigma.alpha[iter],Epsilon[iter],
                    likelihood[iter])
  colnames(save_mcmc) = c("iter",
                          paste("beta",1:p,sep=""),
                          paste("alpha",1:p,sep=""),
                          paste("gamma",1:p,sep=""),
                          "lambda","sigma_gamma","sigma_alpha","epsilon","loglik")

  return(list(post_summary=post_summary, dat = list(X=X,y=y,net=net), save_mcmc = save_mcmc))

}
