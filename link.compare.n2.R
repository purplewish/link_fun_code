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
source('link_fun_code/weights.fun.R')
# ns: sample size, nrep: number of replication; beta0: true parameters;
# min.value, max.value: range of independent variables 
# init: initial values in gev.mle 
# other: initial values of nu in robit model.
# model: the true model assumption, 
#including 'logit','probit','robit','gev','splogit'
# the value of xi, r, nu are specified in model.args 
# inital values of xi r and nu 
### baisc nonlinear form -0.2*(x-3)^2+4
link.compare.n2<- function(model,s0=0,ns,nrep,muv = 0,sdv =1,case=2,model.args = list(),len.newx=200,init.args=list(init=c(0,0),xi0 =1,r0=1,intervalr=c(0.03,10)),spline.control = list(deg = 3,nknots = 10),lamv=seq(1,50,length.out = 20),bound=3,nb=2,iter =1000,weights.arg = 'equal')                                                                                          
{
  ### output ####
  nw <- length(weights.arg)
  
  prmse.ls <- list()
  
  
  mat <- matrix(0,nrep,7)
  colnames(mat)<- c('logit','probit','robit','gev','splogit','pspline(wm)','pspline(m)')
  rmse.mat <- mat
  
  for(w in 1:nw)
  {
    prmse.ls[[w]] <- mat
  }
  
  names(prmse.ls) <- weights.arg
  
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
    
    
    ####new data ##### 
    newdata <- as.matrix(seq(min(x1),max(x1),length.out = len.newx))
    colnames(newdata) <- 'x1'
    
    
    if(case ==1)
    {
      eta0 <- x1^2+x1-3
      eta.new <- newdata^2 + newdata -3
    }
    
    if(case == 2)
    {
      eta0 <- -0.2*(x1-3)^2-0.15*x1+3
      eta.new <- -0.2*(newdata-3)^2-0.15*newdata+3
    }

    
    if(case ==3)
    {
      eta0 <-0.2*(x1+1)^3-2
      eta.new <- 0.2*(newdata+1)^3-2
    }
    
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
    
    prob0[prob0>1]=1
    y0 <- rbinom(ns,size = 1,prob = prob0)
    
    init <- init.args$init
    logit.fit <- glm(y0~x1,family = binomial(link='logit'))
    probit.fit <- glm(y0 ~ x1, family=binomial(link='probit'))
    robit.fit <- robit.pxem(y0 = y0,x0 = x1,beta0 = init,nu0 = init.args$nu0,tol = 1e-3,interval.nu = init.args$interval.nu) 
    gev.fit<- gev.mle.new(y0 = y0,x0 = x1,par0 = c(init,init.args$xi0),maxeval = 50000)
    splogit.fit <- splogit.mle(y0 = y0,x0 = x1,par0 = init,intervalr = init.args$intervalr)
    
    nknots <- spline.control$nknots
    deg <- spline.control$deg
    bs.nc <- nknots+deg-1
    delta0 <- rep(1/bs.nc,bs.nc)
    lam.wm<- pspline.gcv(y0 = y0,x0 = x1,deg = deg, monotone = FALSE,delta0=delta0,nknots = nknots,boundary = range(x1),lamv = lamv,MaxIter=iter)
    
    pspline.fit.wm <- psplinelink(y0 = y0,x0 = x1,deg = deg,lambda = lam.wm, monotone = FALSE,delta0=delta0,nknots = nknots,boundary = range(x1),MaxIter=iter)
    
    lam.m<- pspline.gcv(y0 = y0,x0 = x1,deg = deg, monotone = TRUE,delta0=delta0,nknots = nknots,boundary = range(x1),lamv = lamv,MaxIter=iter)
    
    pspline.fit.m <- psplinelink(y0 = y0,x0 = x1,deg = deg,lambda = lam.m, monotone = TRUE,delta0=delta0,nknots = nknots,boundary = range(x1),MaxIter=iter)
    
    
    
    boundary.n1[j] <-  min(1-gev.fit$est[3]*cbind(1,x1)%*%gev.fit$est[1:2])
    

   
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
    pspline.pred.wm <- predict.pspline(est.obj = pspline.fit.wm,newdata = newdata)
    pspline.pred.m <- predict.pspline(est.obj = pspline.fit.m,newdata = newdata)
    #gam.pred <- predict(gam.fit,newdata = as.data.frame(newdata),type = 'response')
    
 
    
    rmse.mat[j,'logit'] <- sqrt(mean((logit.fit$fitted.values - prob0)^2))
    rmse.mat[j,'probit'] <- sqrt(mean((probit.fit$fitted.values - prob0)^2))
    rmse.mat[j,'robit'] <- sqrt(mean((robit.fit$fitted.values - prob0)^2))
    rmse.mat[j,'gev'] <- sqrt(mean((gev.fit$fitted.values - prob0)^2))
    rmse.mat[j,'splogit'] <- sqrt(mean((splogit.fit$fitted.values - prob0)^2))
    rmse.mat[j,'pspline(wm)']<- sqrt(mean((pspline.fit.wm$fitted.values - prob0)^2))
    rmse.mat[j,'pspline(m)']<- sqrt(mean((pspline.fit.m$fitted.values - prob0)^2))
    
    
    for(w in 1:nw )
    {     
      prmse.ls[[w]][j,'logit'] <- sqrt(sum(weights[,w]*(logit.pred- prob.new)^2))
      prmse.ls[[w]][j,'probit'] <- sqrt(sum(weights[,w]*(probit.pred - prob.new)^2))
      prmse.ls[[w]][j,'robit'] <- sqrt(sum(weights[,w]*(robit.pred- prob.new)^2))
      prmse.ls[[w]][j,'gev'] <- sqrt(sum(weights[,w]*(gev.pred - prob.new)^2))
      prmse.ls[[w]][j,'splogit'] <- sqrt(sum(weights[,w]*(splogit.pred - prob.new)^2))
      prmse.ls[[w]][j,'pspline(wm)'] <- sqrt(sum(weights[,w]*(pspline.pred.wm - prob.new)^2))
      prmse.ls[[w]][j,'pspline(m)'] <- sqrt(sum(weights[,w]*(pspline.pred.m - prob.new)^2))
      
    } 
    
    print(c(j,lam.wm,lam.m))
  }
  

  
  outls <- list(rmse.mat = rmse.mat,prmse.ls=prmse.ls, 
                gr = grv.n1, boundary = boundary.n1,splogit.rv = splogit.rv.n1)
  return(outls)
}
