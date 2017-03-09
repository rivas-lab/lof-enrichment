piemm <- function(betas, ses, priormu0 = 0, priorvar0 = .3^2, priormu1 = 2, priormu2 = -2, priorvar1 = .2^2, priorvar2 = .2^2,burnin = 50, niter = 100, sdconstr = TRUE, nprior = 1, priorsesq=.2^2, igalpha = 1, igbeta = 1, paramfix = FALSE, mu1fix = 2, mu2fix = -2, var1fix = 1, var2fix = 1){
  set.seed(10000)
  y = cbind(betas,ses)
  #Manuel A. Rivas, 03.02.2017
  library(MASS)
  library(mvtnorm)
  igalpha = igalpha
  igbeta = igbeta
  n = nrow(y)
  #print(dim(y))
  ngroups = 3
  na.i = 0
  # priors
  dir.prior = rep(1,ngroups) # dirichlet parameters for pi
  
  # initialize
  gamma = rep(0,n)
  res.mu = matrix(NA,nrow = niter, ncol = ngroups)
  res.mu[1,] = c(priormu0, priormu1, priormu2)
  res.pi = matrix(NA, nrow = niter, ncol = ngroups)
  res.pi[1,] = dir.prior/3
  res.sigma = matrix(NA, nrow = niter, ncol = ngroups)
  res.sigma[1,] = c(priorvar0, priorvar1, priorvar2)
  
  # probability of variant membership for each allele
  m = matrix(0, nrow = n, ncol = ngroups)
  
  # allele membership
  res.gamma = matrix(0, ncol = ngroups, nrow = n)
  
  for(iter in 2:niter){
    #update labelings gamma 
    for(i in 1:ngroups){ 
      m[,i] = log(res.pi[iter-1,i]) + dnorm(betas,res.mu[iter-1,i],sd = sqrt(res.sigma[iter-1,i] + ses^2),log = TRUE)
    }
    m = exp(m - apply(m,1,max))
    m = m/apply(m,1,sum)
    nanum = 0
    for(i in 1:n){
      if(is.na(m[i,1])){ m[i,] = res.pi[1,];
      m[i,] = m[i,]/sum(m[i,])
      nanum = nanum + 1
      }
      sample(1:ngroups, size = 1, prob = m[i,])-> gamma[i]}
    #if(iter%%100==0) print(paste(iter," iterations"))
    if(iter >= burnin){
      for(i in 1:ngroups){
        res.gamma[which(gamma==i),i]=1+res.gamma[which(gamma==i),i]
      }
    }
    # update pi from dirichlet
    for(i in 1:ngroups) {
      pi[i] = rgamma(1,shape=dir.prior[i] + sum(gamma==i), scale = 1)
    }
    pi = pi/sum(pi)
    res.pi[iter,] = pi
    
    # update sigma
    for(i in 1:ngroups){
      which(gamma==i) -> ind
      if(length(ind)>0){
        S = sum((betas[ind] - res.mu[iter-1,i])^2)
        if(i == 1){
          res.sigma[iter,i] = priorvar0
          #res.sigma[iter,i] = 1/(rgamma(1,igalpha + 0.5*length(ind),scale = 1/(igbeta + 0.5*S)))
          
        }
        else{
          if(i == 3){
            if(sdconstr){
              res.sigma[iter, 3] = priorvar2
#              res.sigma[iter,i] = res.sigma[iter,2]
            }else{
              res.sigma[iter,i] = 1/(rgamma(1,shape = igalpha + 0.5*length(ind),scale = 1/(igbeta + 0.5*S)))
            }
          }
          else{
            res.sigma[iter,i] = 1/(rgamma(1,shape = igalpha + 0.5*length(ind),scale = 1/(igbeta + 0.5*S)))
          }
        }
      }
      if(paramfix){
        res.sigma[iter,2] = var1fix
        res.sigma[iter,3] = var2fix
      }
    }
    
    # update mu
    for(i in 1:ngroups){
      which(gamma==i)->ind
      if(length(ind)>0){
        mu.prior = res.mu[1,i]
        sigma.prior = res.sigma[1,i]
        a = (mu.prior/sigma.prior + sum(betas[ind])/res.sigma[iter,i])/(1/sigma.prior + length(ind)/res.sigma[iter,i])
        b = 1/(1/sigma.prior + length(ind)/res.sigma[iter,i])
        if(i == 1){
          res.mu[iter,i] = priormu0
        }else{
        # if((i == 3) & (sdconstr == TRUE) ){
        #   res.mu[iter,i] = priormu2
        # } 
        #  else{
          mu.star = rnorm(1,a,sd=sqrt(b))
          if(is.na(mu.star) | is.na(a)){
            mu.star = rnorm(1,0,1)
            na.i = na.i + 1
          }
          if((mu.star < 0) & (i==2)){
            res.mu[iter,i] = abs(mu.star)
          }else{
            res.mu[iter,i] = mu.star
          }
          
        }
      }
      if(paramfix){
        #   if(paramfix){
        #      res.mu[iter,2] = mu1.fix
        #      res.mu[iter,3] = mu2.fix
        #    }
      }
    }
  }
  l95pi = rep(0,ngroups)
  u95pi = rep(0,ngroups)
  meanpi = rep(0,ngroups)
  l95mu = rep(0,ngroups)
  u95mu = rep(0,ngroups)
  meanmu = rep(0,ngroups)
  l95sd = rep(0,ngroups)
  meansd = rep(0,ngroups)
  u95sd = rep(0,ngroups)
  for(i in 1:ngroups){
    s.pi =  sort(res.pi[(burnin+1):niter,i])
    s.mu = sort(res.mu[(burnin+1):niter,i])
    s.var = sort(res.sigma[(burnin+1):niter,i])
    l95pi[i] = quantile(s.pi, .025)
    u95pi[i] = quantile(s.pi, .975)
    l95mu[i] = quantile(s.mu, .025)
    u95mu[i] = quantile(s.mu, .975)
    l95sd[i] = sqrt(quantile(s.var, .025))
    u95sd[i] = sqrt(quantile(s.var, .975))
    meansd[i] = sqrt(mean(s.var))
    meanpi[i] = mean(s.pi)
    meanmu[i] = mean(s.mu)
  }
  # print(na.i)
  # print("==========================SUMMARY==========================")
  # print(paste(c("proportion: ", meanpi), sep = ' '))
  # print(paste(c("l95 proportion: ", l95pi), sep = ' '))
  # print(paste(c("u95 proportion: ", u95pi), sep = ' '))
  # print(paste(c("mu: ", meanmu), sep = ' '))
  # print(paste(c("l95 mu: ", l95mu), sep = ' '))
  # print(paste(c("u95 mu: ", u95mu), sep = ' '))
  # print(paste(c("sd: ", meansd), sep = ' '))
  # print(paste(c("l95 sd: ", l95sd), sep = ' '))
  # print(paste(c("u95 sd: ", u95sd), sep = ' '))
  x = list(res.mu, res.pi, res.sigma, res.gamma,meansd,meanpi,meanmu,l95sd,l95pi,l95mu,u95sd,u95pi,u95mu)
  return(x)
}
