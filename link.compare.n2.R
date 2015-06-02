### link comparison based on two covariates, nonlinear model ####
## for GEV use uniform distribution####
library(evd)
library(nloptr)
library(bgeva)
library(mgcv)
source('link_fun_code/robit.em.R')
source('link_fun_code/splogit.mle.R')
source('link_fun_code/psplinelink.R')
source('link_fun_code/gev.mle.R')
# ns: sample size, nrep: number of replication; beta0: true parameters;
# min.value, max.value: range of independent variables 
# init: initial values in gev.mle 
# other: initial values of nu in robit model.
# model: the true model assumption, 
#including 'logit','probit','robit','gev','splogit'
# the value of xi, r, nu are specified in model.args 
# inital values of xi r and nu 
### baisc nonlinear form -0.2*(x-3)^2+4
link.compare.n2<- function(model,s0=0,ns,nrep,muv = 0,sdv =1,model.args = list(),len.newx=200,init.args=list(init=c(0,0),xi0 =1,r0=1,intervalr=c(0.03,10)),spline.control = list(deg = 3,nknots = 10),lamv=seq(1,50,length.out = 20),bound=3,nb=2,iter =1000)                                                                                                                                                                                                                                                                                                  
{
  ### output ####
  rmse.logit <- rmse.probit <- rmse.gev <- rmse.robit <- rmse.splogit <- rmse.pspline  <- rep(0,nrep)
  prmse.logit <- prmse.probit <- prmse.gev <- prmse.robit <- prmse.splogit <- prmse.pspline <- rep(0,nrep)
  
  
  #### gradient ####
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
    x1 <- rnorm(4*ns,muv,sd = sdv)
    x1 <- (x1[x1 >= lowp & x1 <= upp])[1:ns]
    eta0 <- -0.2*(x1-3)^2+4
    
    ####new data ##### 
    newdata <- as.matrix(seq(min(x1),max(x1),length.out = len.newx))
    colnames(newdata) <- 'x1'
    eta.new <- -0.2*(newdata-3)^2+4
    
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
      locv <- model.args$locv
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
    logit.fit <- glm(y0~x1,family = binomial(link='logit'))
    probit.fit <- glm(y0 ~ x1, family=binomial(link='probit'))
    robit.fit <- robit.pxem(y0 = y0,x0 = x1,beta0 = init,nu0 = init.args$nu0,tol = 1e-3) 
    gev.fit<- gev.mle.new(y0 = y0,x0 = x1,par0 = c(init,init.args$xi0),maxeval = 50000)
    splogit.fit <- splogit.mle(y0 = y0,x0 = x1,par0 = init,intervalr = init.args$intervalr)
    
    nknots <- spline.control$nknots
    deg <- spline.control$deg
    bs.nc <- nknots+deg-1
    delta0 <- rep(1/bs.nc,bs.nc)
    lam<- pspline.gcv(y0 = y0,x0 = x1,deg = deg, monotone = FALSE,delta0=delta0,nknots = nknots,boundary = range(x1),lamv = lamv,MaxIter=iter)
    
    pspline.fit <- psplinelink(y0 = y0,x0 = x1,deg = deg,lambda = lam, monotone = FALSE,delta0=delta0,nknots = nknots,boundary = range(x1),MaxIter=iter)
    
    
    
    boundary.n1[j] <-  min(1-gev.fit$est[3]*cbind(1,x1)%*%gev.fit$est[1:2])
    
    rmse.logit[j] <- sqrt(mean((logit.fit$fitted.values - prob0)^2))
    rmse.probit[j] <- sqrt(mean((probit.fit$fitted.values - prob0)^2))
    rmse.robit[j] <- sqrt(mean((robit.fit$fitted.values - prob0)^2))
    rmse.gev[j] <- sqrt(mean((gev.fit$fitted.values - prob0)^2))
    rmse.splogit[j] <- sqrt(mean((splogit.fit$fitted.values - prob0)^2))
    rmse.pspline[j] <- sqrt(mean((pspline.fit$fitted.values - prob0)^2))
    #rmse.gam[j] <- sqrt(mean((gam.fit$fitted.values - prob0)^2))
    
    
    grv.n1[j,] <- gev.fit$gr
    
    splogit.rv.n1[j,] <- splogit.fit$gr
    
    
    
    #### ----------------------------predict for new data ------------------------------####
    logit.pred <- predict(logit.fit,newdata = as.data.frame(newdata),type = 'response')
    probit.pred <- predict(probit.fit,newdata = as.data.frame(newdata),type = 'response')
    robit.pred <- predict.pxem(est.obj = robit.fit,newdata = newdata)
    gev.pred <- predict.gev.new(est.obj = gev.fit,newdata = newdata)
    splogit.pred <- predict.splogit(est.obj = splogit.fit,newdata = newdata)
    pspline.pred <- predict.pspline(est.obj = pspline.fit,newdata = newdata)
    #gam.pred <- predict(gam.fit,newdata = as.data.frame(newdata),type = 'response')
    
    prmse.logit[j] <- sqrt(mean((logit.pred- prob.new)^2))
    prmse.probit[j] <- sqrt(mean((probit.pred - prob.new)^2))
    prmse.robit[j] <- sqrt(mean((robit.pred- prob.new)^2))
    prmse.gev[j] <- sqrt(mean((gev.pred - prob.new)^2))
    prmse.splogit[j] <- sqrt(mean((splogit.pred - prob.new)^2))
    prmse.pspline[j] <- sqrt(mean((pspline.pred - prob.new)^2))
    #prmse.gam[j] <- sqrt(mean((as.numeric(gam.pred) - prob.new)^2))
    
    
    print(j)
  }
  
  rmse.mat <- cbind(rmse.logit,rmse.probit,rmse.robit,
                    rmse.gev,rmse.splogit,rmse.pspline)
  
  prmse.mat <- cbind(prmse.logit,prmse.probit,prmse.robit,
                     prmse.gev,prmse.splogit,prmse.pspline)
  
  outls <- list(rmse.mat = rmse.mat,prmse.mat=prmse.mat, 
                gr = grv.n1, boundary = boundary.n1,splogit.rv = splogit.rv.n1)
  return(outls)
}
