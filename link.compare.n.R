### link comparison based on two covariates, nonlinear model ####
## for GEV use uniform distribution####
library(evd)
library(nloptr)
library(bgeva)
library(mgcv)
source('link_fun_code/robit.em.R')
source('link_fun_code/splogit.mle.R')
source('link_fun_code/psplinelink1.R')
source('link_fun_code/gev.mle.R')
# ns: sample size, nrep: number of replication; beta0: true parameters;
# min.value, max.value: range of independent variables 
# init: initial values in gev.mle 
# other: initial values of nu in robit model.
# model: the true model assumption, 
#including 'logit','probit','robit','gev','splogit'
# the value of xi, r, nu are specified in model.args 
# inital values of xi r and nu 
link.compare.n<- function(model,s0=0,ns,nrep,min.value,max.value,muv = c(0,1),sdv =c(sqrt(3),sqrt(3)),  model.args = list(),init.args=list(init=c(0,0,0),xi0 =1,r0=1,intervalr=c(0.03,10)),
spline.control = list(deg = 3,nknots = 10,beta0s=c(1,1)),lamv=seq(5,100,length.out = 20),bound=2.5)                                                                                                                                                                                                                                                                                                  
{
  ### output ####
  mse.logit <- mse.probit <- mse.gev <- mse.robit <- mse.splogit <- mse.pspline <-  mse.gam <- rep(0,nrep)
  
  #### aic ##### 
  #aic.logit <- aic.probit <- aic.gev <- aic.gev.new <- aic.robit <- aic.splogit <- aic.pspline  <- rep(0,nrep)
  
  #### gradient ####
  nb <- 3
  grv.n1 <- matrix(0,nrep,nb+1)
  splogit.rv.n1 <- matrix(0,nrep,nb)
  boundary.n1 <- rep(0,nrep)
  
  deg <- spline.control$deg
  nknots <- spline.control$nknots
  beta0s <- spline.control$beta0s
  
  for(s in 1:nrep)
  {
    set.seed(s+s0)
    x1 <- runif(ns,min = -4,max = -1)
    x2 <- runif(ns,min = 2 ,max = 4)
    yita0 <- (x1+x2)^2+ x1
#     x1 <- rnorm(4*ns,muv[1],sd = sdv[1])
#     x2 <- rnorm(4*ns,muv[2], sd = sdv[2])
#     x1 <- x1[abs(x1) <= bound][1:ns]
#     x2 <- x2[abs(x2) <= bound][1:ns]
#     yita0 <- -1.5+0.5*(x1+x2)^2+x1*x2
    
    if(model == 'logit')
    {
      prob0 <- exp(yita0)/(1+exp(yita0))
    }
    
    if(model == 'probit')
    {
      prob0 <- pnorm(q = yita0)
    }
    
    if(model=='robit')
    {
      nu <- model.args$nu
      prob0 <- pt(q = yita0,df = nu)
    }
    if(model=='gev')
    {
      xi <- model.args$xi
#     x1 <- runif(ns,min = min.value,max = max.value)
#     x2 <- rbinom(ns,size=1,prob=0.5)
#     yita0 <- cbind(1,x1,x2)%*%betav
      prob0 <- 1-pgev(-yita0,loc = 0,scale = 1,shape = xi)
    }
    
    if(model == 'splogit')
    {
#       x1 <- rnorm(2*ns,0,sd = sdv)
#       x1 <- (x1[abs(x1) <= bound])[1:ns]
#       x2 <- rbinom(ns,size=1,prob=0.5)

      
      splogit.link <- function(yita0,r0)
      {
        if(r0 >0 & r0 <=1)
        {
          prob0 <- exp(yita0)/((1+exp(yita0/r0))^r0)
        }
        
        if(r0 > 1)
        {
          yita.new <- -r0*yita0
          prob0 <- 1-(exp(yita.new)/(1+exp(yita.new)))^(1/r0)
        }
        return(prob0)
      }
      
      r <- model.args$r
      prob0 <- splogit.link(yita0,r)
      
    }
    
    y0 <- rbinom(ns,size = 1,prob = prob0)
    xmat0 <- cbind(x1,x2)
    
    init <- init.args$init
    logit.fit <- glm(y0~x1+x2,family = binomial(link='logit'))
    probit.fit <- glm(y0 ~ x1+x2, family=binomial(link='probit'))
    robit.fit <- robit.pxem(y0 = y0,x0 = xmat0,beta0 = init,nu0 = init.args$nu0,tol = 1e-3) 
    gev.fit<- gev.mle.new(y0 = y0,x0 = xmat0,par0 = c(init,init.args$xi0),maxeval = 50000)
    splogit.fit <- splogit.mle(y0 = y0,x0 = xmat0,par0 = init,intervalr = init.args$intervalr)
    
    
    lam<- pspline.gcv(y0 = y0,xmat = xmat0,qv=1,monotone = TRUE,nknots = 10,beta0 = c(1,1),lamv = lamv,d.value = 9,MaxIter = 1000)
    
    pspline.fit <- psplinelink1(y0 = y0,xmat = xmat0,qv=1,monotone = TRUE,nknots = 10,beta0 = beta0s,lambda=lam,d.value = 9,MaxIter = 1000)
    
    gam.fit<- gam(y0~s(x1)+x2,family = binomial(link = 'logit'))
    
    boundary.n1[s] <-  min(1-gev.fit$est[4]*cbind(1,xmat0)%*%gev.fit$est[1:3])
    
    mse.logit[s] <- mean((logit.fit$fitted.values - prob0)^2)
    mse.probit[s] <- mean((probit.fit$fitted.values - prob0)^2)
    mse.robit[s] <- mean((robit.fit$fitted.values - prob0)^2)
    mse.gev[s] <- mean((gev.fit$fitted.values - prob0)^2)
    mse.splogit[s] <- mean((splogit.fit$fitted.values - prob0)^2)
    mse.pspline[s] <- mean((pspline.fit$fitted.values - prob0)^2)
    mse.gam[s] <- mean((gam.fit$fitted.values - prob0)^2)
    
    grv.n1[s,] <- gev.fit$gr
    
    splogit.rv.n1[s,] <- splogit.fit$gr
    
    print(s)
  }
  
  mse.mat <- cbind(mse.logit,mse.probit,mse.robit,
                   mse.gev,mse.splogit,mse.gam,mse.pspline)
  
  outls <- list(mse.mat = mse.mat, 
                gr = grv.n1, boundary = boundary.n1,splogit.rv = splogit.rv.n1)
  return(outls)
}


