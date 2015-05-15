### link comparison based on two covariates, one is N(0,3), another one is binary####
## for GEV use uniform distribution####
library(evd)
library(nloptr)
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
## bound*sd 
link.compare.b<- function(model,s0=0,ns,nrep,muv=0,sdv =1,bound=3,len.newx=200,
                        model.args=list(beta0=c(0,1,1)),
                        init.args=list(init=c(0,0,0),xi0 =1,r0=1,nu0=1,intervalr=c(0.03,10)),
                        spline.control = list(deg = 3,nknots = 10),lamv=seq(1,20,length.out = 10),iter=50)                                                                                                                                                                                                                                                                                                  
{
  
  ### output ####
  rmse.logit <- rmse.probit <- rmse.gev <- rmse.robit <- rmse.splogit <- rmse.pspline <-  rmse.gam <- rep(0,nrep)
  prmse.logit <- prmse.probit <- prmse.gev <- prmse.robit <- prmse.splogit <- prmse.pspline <-  prmse.gam <- rep(0,nrep)
  
  rrmse.logit <- rrmse.probit <- rrmse.gev <- rrmse.robit <- rrmse.splogit <- rrmse.pspline <-  rrmse.gam <- rep(0,nrep)
  prrmse.logit <- prrmse.probit <- prrmse.gev <- prrmse.robit <- prrmse.splogit <- prrmse.pspline <-  prrmse.gam <- rep(0,nrep)
  
  #### aic ##### 
  #aic.logit <- aic.probit <- aic.gev <- aic.gev.new <- aic.robit <- aic.splogit <- aic.pspline  <- rep(0,nrep)
  

  betav <- model.args$beta0
  nb <- length(betav)
  ## gradient ###
 
  grv.n1 <- matrix(0,nrep,nb+1)
  splogit.rv.n1 <- matrix(0,nrep,nb)
  boundary.n1 <- rep(0,nrep)
  
  deg <- spline.control$deg
  nknots <- spline.control$nknots
  
  lowp <- muv - bound*sdv
  upp <- muv + bound*sdv
  
  for(j in 1:nrep)
  {
    set.seed(j+s0)
    x1 <- rnorm(2*ns,muv,sd = sdv)
    x1 <- (x1[x1 >= lowp & x1 <= upp])[1:ns]
    x2 <- rbinom(ns,size=1,prob=0.5)
    eta0 <- cbind(1,x1,x2)%*%betav
    xmat0 <- cbind(x1,x2)
    
    ### new data for prediction #####
    newdata <- as.matrix(expand.grid(seq(min(x1),max(x2),length.out = len.newx),c(0,1)))
    colnames(newdata) <- colnames(xmat0)
    eta.new <- cbind(1,newdata)%*%betav
    
    if(model == 'logit')
    {
      prob0 <- exp(eta0)/(1+exp(eta0))
      prob.new <- exp(eta.new)/(1+exp(eta.new))
    }
    
    if(model == 'probit')
    {
      prob0 <- pnorm(q = eta0)
      prob.new <- pnorm(q=eta.new)
    }
    
    if(model=='robit')
    {
      nu <- model.args$nu
      prob0 <- pt(q = eta0,df = nu)
      prob.new <- pt(q = eta.new,df = nu)
    }
    if(model=='gev')
    {
      xi <- model.args$xi
      locv <- model.args$loc
      prob0 <- 1-pgev(-eta0,loc = locv,scale = 1,shape = xi)
      prob.new <- 1-pgev(-eta.new,loc = 0,scale = 1,shape = xi)
    }
    
    if(model == 'splogit')
    {
      
      splogit.link <- function(eta0,r0)
      {
        if(r0 >0 & r0 <=1)
        {
          prob0 <- exp(eta0)/((1+exp(eta0/r0))^r0)
        }
        
        if(r0 > 1)
        {
          eta.new <- -r0*eta0
          prob0 <- 1-(exp(eta.new)/(1+exp(eta.new)))^(1/r0)
        }
        return(prob0)
      }
      
      r <- model.args$r
      prob0 <- splogit.link(eta0,r)
      prob.new <- splogit.link(eta.new,r)
      
    }
    
    y0 <- rbinom(ns,size = 1,prob = prob0)

    
    init <- init.args$init
    logit.fit <- glm(y0~x1+x2,family = binomial(link='logit'))
    probit.fit <- glm(y0 ~ x1+x2, family=binomial(link='probit'))
    robit.fit <- robit.pxem(y0 = y0,x0 = xmat0,beta0 = init,nu0 = init.args$nu0,tol = 1e-3) 
    gev.fit<- gev.mle.new(y0 = y0,x0 = xmat0,par0 = c(init,init.args$xi0),maxeval = 50000)
    splogit.fit <- splogit.mle(y0 = y0,x0 = xmat0,par0 = init,intervalr = init.args$intervalr)
    
  
    lam<- pspline.gcv(y0 = y0,xmat = xmat0,qv=1,catv = 'x2',monotone = TRUE,nknots = nknots,beta0 = c(1,1),MaxIter = 500,lamv = lamv)
    
    pspline.fit <- psplinelink1(y0 = y0,xmat = xmat0,qv=1,catv = 'x2',monotone = TRUE,nknots = nknots,beta0 = c(1,1),lambda=20,MaxIter = iter)
    
    gam.fit<- gam(y0~s(x1)+x2,family = binomial(link = 'logit'))
    
    boundary.n1[j] <-  min(1-gev.fit$est[4]*cbind(1,xmat0)%*%gev.fit$est[1:3])
  
    rmse.logit[j] <- sqrt(mean((logit.fit$fitted.values - prob0)^2))
    rmse.probit[j] <- sqrt(mean((probit.fit$fitted.values - prob0)^2))
    rmse.robit[j] <- sqrt(mean((robit.fit$fitted.values - prob0)^2))
    rmse.gev[j] <- sqrt(mean((gev.fit$fitted.values - prob0)^2))
    rmse.splogit[j] <- sqrt(mean((splogit.fit$fitted.values - prob0)^2))
    rmse.pspline[j] <- sqrt(mean((pspline.fit$fitted.values - prob0)^2))
    rmse.gam[j] <- sqrt(mean((gam.fit$fitted.values - prob0)^2))
    
    rrmse.logit[j] <- sqrt(mean((logit.fit$fitted.values - prob0)^2/prob0^2))
    rrmse.probit[j] <- sqrt(mean((probit.fit$fitted.values - prob0)^2/prob0^2))
    rrmse.robit[j] <- sqrt(mean((robit.fit$fitted.values - prob0)^2/prob0^2))
    rrmse.gev[j] <- sqrt(mean((gev.fit$fitted.values - prob0)^2/prob0^2))
    rrmse.splogit[j] <- sqrt(mean((splogit.fit$fitted.values - prob0)^2/prob0^2))
    rrmse.pspline[j] <- sqrt(mean((pspline.fit$fitted.values - prob0)^2/prob0^2))
    rrmse.gam[j] <- sqrt(mean((gam.fit$fitted.values - prob0)^2/prob0^2))
    
     
    grv.n1[j,] <- gev.fit$gr

    splogit.rv.n1[j,] <- splogit.fit$gr
    
    
    
    #### ----------------------------predict for new data ------------------------------####
    logit.pred <- predict(logit.fit,newdata = as.data.frame(newdata),type = 'response')
    probit.pred <- predict(probit.fit,newdata = as.data.frame(newdata),type = 'response')
    robit.pred <- predict.pxem(est.obj = robit.fit,newdata = newdata)
    gev.pred <- predict.gev.new(est.obj = gev.fit,newdata = newdata)
    splogit.pred <- predict.splogit(est.obj = splogit.fit,newdata = newdata)
    pspline.pred <- predict.pspline1(est.obj = pspline.fit,newdata = newdata)
    gam.pred <- predict(gam.fit,newdata = as.data.frame(newdata),type = 'response')
    
    prmse.logit[j] <- sqrt(mean((logit.pred- prob.new)^2))
    prmse.probit[j] <- sqrt(mean((probit.pred - prob.new)^2))
    prmse.robit[j] <- sqrt(mean((robit.pred- prob.new)^2))
    prmse.gev[j] <- sqrt(mean((gev.pred - prob.new)^2))
    prmse.splogit[j] <- sqrt(mean((splogit.pred - prob.new)^2))
    prmse.pspline[j] <- sqrt(mean((pspline.pred - prob.new)^2))
    prmse.gam[j] <- sqrt(mean((as.numeric(gam.pred) - prob.new)^2))
    
    prrmse.logit[j] <- sqrt(mean((logit.pred- prob.new)^2/prob.new^2))
    prrmse.probit[j] <- sqrt(mean((probit.pred - prob.new)^2/prob.new^2))
    prrmse.robit[j] <- sqrt(mean((robit.pred- prob.new)^2/prob.new^2))
    prrmse.gev[j] <- sqrt(mean((gev.pred - prob.new)^2/prob.new^2))
    prrmse.splogit[j] <- sqrt(mean((splogit.pred - prob.new)^2/prob.new^2))
    prrmse.pspline[j] <- sqrt(mean((pspline.pred - prob.new)^2/prob.new^2))
    prrmse.gam[j] <- sqrt(mean((as.numeric(gam.pred) - prob.new)^2/prob.new^2))
    
    
    print(j)
  }
  
  rmse.mat <- cbind(rmse.logit,rmse.probit,rmse.robit,
                   rmse.gev,rmse.splogit,rmse.gam,rmse.pspline)
  
  rrmse.mat <- cbind(rrmse.logit,rrmse.probit,rrmse.robit,
                    rrmse.gev,rrmse.splogit,rrmse.gam,rrmse.pspline)
  
  prmse.mat <- cbind(prmse.logit,prmse.probit,prmse.robit,
                    prmse.gev,prmse.splogit,prmse.gam,prmse.pspline)
  
  prrmse.mat <- cbind(prrmse.logit,prrmse.probit,prrmse.robit,
                     prrmse.gev,prrmse.splogit,prrmse.gam,prrmse.pspline)
  
  outls <- list(rmse.mat = rmse.mat, rrmse.mat = rrmse.mat,
                prmse.mat=prmse.mat,prrmse.mat = prrmse.mat,
                gr = grv.n1, boundary = boundary.n1,splogit.rv = splogit.rv.n1)
  return(outls)
}



