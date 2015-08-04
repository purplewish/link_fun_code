### link comparison based on two covariates, one is N(0,3), another one is binary####
## for GEV use uniform distribution####
library(evd)
library(nloptr)
library(mgcv)
source('link_fun_code/robit.em.R')
source('link_fun_code/splogit.mle.R')
source('link_fun_code/psplinelink3.R')
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


<<<<<<< HEAD
link.compare.b3<- function(model,s0=0,ns,nrep,muv=0,sdv =1,bound=3,len.newx=200,weights.arg='equal',model.args=list(beta0=c(0,1,1)),init.args=list(init=c(0,0,0),xi0 =1,r0=1,nu0=1,intervalr=c(0.03,10)), spline.control = list(deg = 3,nknots = 10,dd=1),lamv=seq(5,50,length.out = 20),iter=2000)                             
=======
link.compare.b3<- function(model,s0=0,ns,nrep,muv=0,sdv =1,bound=3,len.newx=200,weights.arg='equal',model.args=list(beta0=c(0,1,1)),init.args=list(init=c(0,0,0),xi0 =1,r0=1,nu0=1,intervalr=c(0.03,10),interval.nu=c(0.1,10)), spline.control = list(deg = 3,nknots = 10,dd=1),lamv=seq(5,50,length.out = 20),iter=2000)                                                                                                                                                                                                  
>>>>>>> 05f78c7f19cb7c6ad6fe4be334609ee509aa2390
{
  ### output ####
  nw <- length(weights.arg)

  prmse.ls <- list()
  
  mat <- matrix(0,nrep,7)
  colnames(mat)<- c('logit','probit','robit','gev','splogit','pspline','gam')
  rmse.mat <- mat
  
  for(w in 1:nw)
  {
   prmse.ls[[w]] <- mat
  }
  
   names(prmse.ls) <- weights.arg
  
  
  betav <- model.args$beta0
  nb <- length(betav)
  ## gradient ###
  
  grv.n1 <- matrix(0,nrep,nb+1)
  splogit.rv.n1 <- matrix(0,nrep,nb)
  boundary.n1 <- rep(0,nrep)
  

  
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
    
    
    dd <- spline.control$dd
    init <- init.args$init
    
    
    deg <- spline.control$deg
    nknots <- spline.control$nknots
    
    logit.fit <- glm(y0~x1+x2,family = binomial(link='logit'))
    probit.fit <- glm(y0 ~ x1+x2, family=binomial(link='probit'))
    robit.fit <- robit.pxem(y0 = y0,x0 = xmat0,beta0 = init,nu0 = init.args$nu0,tol = 1e-3,interval.nu=init.args$interval.nu) 
    gev.fit<- gev.mle.new(y0 = y0,x0 = xmat0,par0 = c(init,init.args$xi0),maxeval = 50000)
    splogit.fit <- splogit.mle(y0 = y0,x0 = xmat0,par0 = init,intervalr = init.args$intervalr)
    
    
    lam<- pspline.gcv3(y0 = y0,xmat = xmat0,qv=1,catv = 'x2',monotone = TRUE,nknots = nknots,beta0 = c(1,1),MaxIter = iter,lamv = lamv,dd=dd)
    
    pspline.fit <- psplinelink3(y0 = y0,xmat = xmat0,qv=1,catv = 'x2',monotone = TRUE,nknots = nknots,beta0 = c(1,1),lambda=lam,MaxIter = iter,dd=dd)
    
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
 
    #### ----------------------------predict for new data ------------------------------####
    logit.pred <- predict(logit.fit,newdata = as.data.frame(newdata),type = 'response')
    probit.pred <- predict(probit.fit,newdata = as.data.frame(newdata),type = 'response')
    robit.pred <- predict.pxem(est.obj = robit.fit,newdata = newdata)
    gev.pred <- predict.gev.new(est.obj = gev.fit,newdata = newdata)
    splogit.pred <- predict.splogit(est.obj = splogit.fit,newdata = newdata)
    pspline.pred <- predict.pspline3(est.obj = pspline.fit,newdata = newdata)
    gam.pred <- predict(gam.fit,newdata = as.data.frame(newdata),type = 'response')
    
    
    rmse.mat[j,'logit'] <- sqrt(mean((logit.fit$fitted.values - prob0)^2))
    rmse.mat[j,'probit'] <- sqrt(mean((probit.fit$fitted.values - prob0)^2))
    rmse.mat[j,'robit'] <- sqrt(mean((robit.fit$fitted.values - prob0)^2))
    rmse.mat[j,'gev'] <- sqrt(mean((gev.fit$fitted.values - prob0)^2))
    rmse.mat[j,'splogit'] <- sqrt(mean((splogit.fit$fitted.values - prob0)^2))
    rmse.mat[j,'pspline']<- sqrt(mean((pspline.fit$fitted.values - prob0)^2))
    rmse.mat[j,'gam'] <- sqrt(mean((gam.fit$fitted.values - prob0)^2))
    
      
    for(w in 1:nw )
    {     
      prmse.ls[[w]][j,'logit'] <- sqrt(sum(weights[,w]*(logit.pred- prob.new)^2))
      prmse.ls[[w]][j,'probit'] <- sqrt(sum(weights[,w]*(probit.pred - prob.new)^2))
      prmse.ls[[w]][j,'robit'] <- sqrt(sum(weights[,w]*(robit.pred- prob.new)^2))
      prmse.ls[[w]][j,'gev'] <- sqrt(sum(weights[,w]*(gev.pred - prob.new)^2))
      prmse.ls[[w]][j,'splogit'] <- sqrt(sum(weights[,w]*(splogit.pred - prob.new)^2))
      prmse.ls[[w]][j,'pspline'] <- sqrt(sum(weights[,w]*(pspline.pred - prob.new)^2))
      prmse.ls[[w]][j,'gam'] <- sqrt(sum(weights[,w]*(as.numeric(gam.pred) - prob.new)^2))
      
    } 
    print(c(j,lam))
  }
  


   
  outls <- list(rmse.mat = rmse.mat, 
                prmse.ls=prmse.ls,
                gr = grv.n1, boundary = boundary.n1,splogit.rv = splogit.rv.n1)
  return(outls)
}


