### link comparison based on two covariates, one is N(0,3), another one is binary####
## for GEV use uniform distribution####
library(evd)
library(nloptr)
library(mgcv)
library(quadprog)
source('link_fun_code/robit.em.R')
source('link_fun_code/splogit.mle.R')
source('link_fun_code/psplinelink5.R')
source('link_fun_code/gev.mle.R')
source('link_fun_code/weights.fun.R')
# ns: sample size, nrep: number of replication; beta0: true parameters;
# min.value, max.value: range of independent variables 
# init: initial values in gev.mle 
# other: initial values of nu in robit model.
# model: the true model assumption, 
#including 'logit','probit','robit','gev','splogit'
# the value of xi, r, nu are specified in model.args 
# inital values of xi r and nu 
## bound*sd 
### weights.arg is a vector to specify what kinds of weights are given in the prediction. it can be "equal",'both','left','right'

ntruncated <- 3
truncate.fun<- function(x,cutoff.value)
{
  xt <- unlist(lapply(unlist(lapply(x,function(x) max(x,cutoff.value[1]))),function(x) min(x,cutoff.value[2])))
  return(xt)
}
link.compare.b5<- function(model,s0=0,ns,nrep,muv=-0.5,sdv =1,bound=3,len.newx=200,weights.arg='equal',model.args=list(beta0=c(0,1,1)),init.args=list(init=c(0,0,0),xi0 =1,r0=1,nu0=1,intervalr=c(0.03,10)), spline.control = list(deg = 3,nknots = 11,qv=0.95),lamv=seq(5,50,length.out = 20),truncated.arg = list(lower=0.25,upper=0.75),iter=2000)                             
{
  ### output ####
  nw <- length(weights.arg)

  prmse.ls <- list()
  truncated.ls <- list()
  
  mat <- matrix(0,nrep,8)
  colnames(mat)<- c('logit','probit','robit','gev','splogit','pspline(wm)','pspline(m)','gam')
  rmse.mat <- mat
  
  for(w in 1:nw)
  {
    prmse.ls[[w]] <- mat
  }
  
  for(w in 1:3)
  {
    truncated.ls[[w]] <- mat
  }
  
  names(prmse.ls) <- weights.arg
  names(truncated.ls) <- c("greater","less","both")
  
  betav <- model.args$beta0
  nb <- length(betav)
  ## gradient ###
  
  grv.n1 <- matrix(0,nrep,nb+1)
  splogit.rv.n1 <- matrix(0,nrep,nb)
  boundary.n1 <- rep(0,nrep)
  lam.track <- matrix(0,nrep,2)
  
  
  if(is.null(bound)==FALSE)
  {
    lowp <- muv - bound*sdv
    upp <- muv + bound*sdv
  }
  
  
  for(j in 1:nrep)
  {
    set.seed(j+s0)
    if(is.null(bound)==FALSE)
    {
      x1 <- rnorm(2*ns,muv,sd = sdv)
      x1 <- (x1[x1 >= lowp & x1 <= upp])[1:ns]
    }
    if(is.null(bound)==TRUE)
    {
      x1 <- rnorm(ns,muv,sd=sdv)
    }
    
    x2 <- rbinom(ns,size=1,prob=0.5)
    eta0 <- cbind(1,x1,x2)%*%betav
    xmat0 <- cbind(x1,x2)
    
    ### new data for prediction #####
    newdata <- as.matrix(expand.grid(seq(min(x1),max(x1),length.out = len.newx),c(0,1)))
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
      prob.new <- 1-pgev(-eta.new,loc = locv,scale = 1,shape = xi)
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
    
    deg <- spline.control$deg
    nknots <- spline.control$nknots
    
    logit.fit <- glm(y0~x1+x2,family = binomial(link='logit'))
    probit.fit <- glm(y0 ~ x1+x2, family=binomial(link='probit'))
    robit.fit <- robit.pxem(y0 = y0,x0 = xmat0,beta0 = init,nu0 = init.args$nu0,tol = 1e-3,interval.nu=init.args$interval.nu) 
    gev.fit<- gev.mle.new(y0 = y0,x0 = xmat0,par0 = c(init,init.args$xi0),maxeval = 50000)
    splogit.fit <- splogit.mle(y0 = y0,x0 = xmat0,par0 = init,intervalr = init.args$intervalr)
    
    qv <- spline.control$qv
    lam.m<- pspline.gcv5(y0 = y0,xmat = xmat0,qv=qv,catv = 'x2',monotone = TRUE,nknots = nknots,beta0 = c(1,1),MaxIter = iter,lamv = lamv)

    
    pspline.fit.m <- psplinelink5(y0 = y0,xmat = xmat0,qv=qv,catv = 'x2',monotone = TRUE,nknots = nknots,beta0 = c(1,1),lambda=lam.m,MaxIter = iter)
    
    lam.wm<- pspline.gcv5(y0 = y0,xmat = xmat0,qv=qv,catv = 'x2',monotone = FALSE,nknots = nknots,beta0 = c(1,1),MaxIter = iter,lamv = lamv)
    
    
    pspline.fit.wm <- psplinelink5(y0 = y0,xmat = xmat0,qv=qv,catv = 'x2',monotone = FALSE,nknots = nknots,beta0 = c(1,1),lambda=lam.wm,MaxIter = iter)
    
    
    
    gam.fit<- gam(y0~s(x1)+x2,family = binomial(link = 'logit'))
    
    boundary.n1[j] <-  min(1-gev.fit$est[4]*cbind(1,xmat0)%*%gev.fit$est[1:3])
    
    
    
    grv.n1[j,] <- gev.fit$gr
    
    splogit.rv.n1[j,] <- splogit.fit$gr
    
    weights <- matrix(0,nrow(newdata),nw)
    colnames(weights) <- weights.arg
    for(w in 1:nw)
    {
      weights[,w] <- weights.fun(weights.arg[w],data = newdata[,1])
    }
    
    truncated.mat <- matrix(c(truncated.arg$upper,1,0,truncated.arg$lower,truncated.arg$lower,truncated.arg$upper),nrow=3,byrow=TRUE)
    

    
    #### ----------------------------predict for new data ------------------------------####
    logit.pred <- predict(logit.fit,newdata = as.data.frame(newdata),type = 'response')
    probit.pred <- predict(probit.fit,newdata = as.data.frame(newdata),type = 'response')
    robit.pred <- predict.pxem(est.obj = robit.fit,newdata = newdata)
    gev.pred <- predict.gev.new(est.obj = gev.fit,newdata = newdata)
    splogit.pred <- predict.splogit(est.obj = splogit.fit,newdata = newdata)
    pspline.pred.m <- predict.pspline5(est.obj = pspline.fit.m,newdata = newdata)
    pspline.pred.wm <- predict.pspline5(est.obj = pspline.fit.wm,newdata = newdata)
    gam.pred <- predict(gam.fit,newdata = as.data.frame(newdata),type = 'response')
    
    
    rmse.mat[j,'logit'] <- sqrt(mean((logit.fit$fitted.values - prob0)^2))
    rmse.mat[j,'probit'] <- sqrt(mean((probit.fit$fitted.values - prob0)^2))
    rmse.mat[j,'robit'] <- sqrt(mean((robit.fit$fitted.values - prob0)^2))
    rmse.mat[j,'gev'] <- sqrt(mean((gev.fit$fitted.values - prob0)^2))
    rmse.mat[j,'splogit'] <- sqrt(mean((splogit.fit$fitted.values - prob0)^2))
    rmse.mat[j,"pspline(wm)"] <- sqrt(mean((pspline.fit.wm$fitted.values - prob0)^2))
    rmse.mat[j,'pspline(m)']<- sqrt(mean((pspline.fit.m$fitted.values - prob0)^2))
    rmse.mat[j,'gam'] <- sqrt(mean((gam.fit$fitted.values - prob0)^2))
    
    
    for(w in 1:nw )
    {     
      prmse.ls[[w]][j,'logit'] <- sqrt(sum(weights[,w]*(logit.pred- prob.new)^2))
      prmse.ls[[w]][j,'probit'] <- sqrt(sum(weights[,w]*(probit.pred - prob.new)^2))
      prmse.ls[[w]][j,'robit'] <- sqrt(sum(weights[,w]*(robit.pred- prob.new)^2))
      prmse.ls[[w]][j,'gev'] <- sqrt(sum(weights[,w]*(gev.pred - prob.new)^2))
      prmse.ls[[w]][j,'splogit'] <- sqrt(sum(weights[,w]*(splogit.pred - prob.new)^2))
      prmse.ls[[w]][j,'pspline(wm)'] <- sqrt(sum(weights[,w]*(pspline.pred.wm - prob.new)^2))
      prmse.ls[[w]][j,'pspline(m)'] <- sqrt(sum(weights[,w]*(pspline.pred.m - prob.new)^2))
      prmse.ls[[w]][j,'gam'] <- sqrt(sum(weights[,w]*(as.numeric(gam.pred) - prob.new)^2))
      
    } 
    
    for(w in 1:3)
    {
      truncated.prob <- truncate.fun(prob.new,truncated.mat[w,])
      truncated.ls[[w]][j,"logit"] <- mean(abs(truncate.fun(logit.pred,truncated.mat[w,])- truncated.prob))
      truncated.ls[[w]][j,"probit"] <- mean(abs(truncate.fun(probit.pred,truncated.mat[w,])- truncated.prob))
      truncated.ls[[w]][j,"robit"] <- mean(abs(truncate.fun(robit.pred,truncated.mat[w,])- truncated.prob))
      truncated.ls[[w]][j,"gev"] <- mean(abs(truncate.fun(gev.pred,truncated.mat[w,])- truncated.prob))
      truncated.ls[[w]][j,"splogit"] <- mean(abs(truncate.fun(splogit.pred,truncated.mat[w,])- truncated.prob))
      truncated.ls[[w]][j,"pspline(wm)"] <- mean(abs(truncate.fun(pspline.pred.wm,truncated.mat[w,])- truncated.prob))
      truncated.ls[[w]][j,"pspline(m)"] <- mean(abs(truncate.fun(pspline.pred.m,truncated.mat[w,])- truncated.prob))
      truncated.ls[[w]][j,"gam"] <- mean(abs(truncate.fun(gam.pred,truncated.mat[w,])- truncated.prob))
      
    }
   
    lam.track[j,] <- c(lam.m,lam.wm)
    print(j)
  }
  
  outls <- list(rmse.mat = rmse.mat, 
                prmse.ls=prmse.ls,
                truncated.ls = truncated.ls,
                gr = grv.n1, boundary = boundary.n1,splogit.rv = splogit.rv.n1,lam.track=lam.track)
  return(outls)
}


