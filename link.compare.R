library(evd)
library(nloptr)
library(bgeva)
source('link_fun_code/robit.em.R')
source('link_fun_code/splogit.mle.R')
source('link_fun_code/splinelink.R')
source('link_fun_code/gev.mle.R')
# ns: sample size, nrep: number of replication; beta0: true parameters;
# min.value, max.value: range of independent variables 
# init: initial values in gev.mle 
# other: initial values of nu in robit model.
# model: the true model assumption, 
#including 'logit','probit','robit','gev','splogit'
# the value of xi, r, nu are specified in model.args 
# inital values of xi r and nu 
link.compare<- function(model,s0=0,ns,nrep,min.value,max.value,
                        model.args=list(beta0=c(1,1)),quantile =c(0.1,0.9),
                        init.args=list(init=c(0,0),xi0 =1,r0=1,intervalr=c(0.03,10)),
                        spline.control = list(deg = 3,nknots = 20))                                                                                                                                                                                                                                                                                                  
{
  ### output ####
  mse.logit <- mse.probit <- mse.gev <- mse.gev.new <- mse.robit <- mse.splogit <- mse.pspline <- ks.logit <- ks.probit <- ks.gev <- ks.robit <- ks.splogit <- ks.pspline <- ks.gev.new <- rep(0,nrep)
  max.logit <- max.probit <- max.gev <- max.gev.new <- max.robit <- max.splogit <- max.pspline <- rep(0,nrep) 
  qtle.logit <- qtle.probit <- qtle.gev <- qtle.gev.new <- qtle.robit <- qtle.splogit <- qtle.pspline <- matrix(0,nrep,length(quantile))
  
  
  #### aic ##### 
  aic.logit <- aic.probit <- aic.gev <- aic.gev.new <- aic.robit <- aic.splogit <- aic.pspline  <- rep(0,nrep)
  
  #### gradient ####
  grv1.n1<- grv2.n1 <- matrix(0,nrep,3)
  splogit.rv.n1 <- matrix(0,nrep,2)
  boundary1.n1 <- boundary2.n1 <- rep(0,nrep)
  
  deg <- spline.control$deg
  nknots <- spline.control$nknots
  
  betav <- model.args$beta0
  
  for(s in 1:nrep)
  {
    set.seed(s+s0)
    x0 <- sort(runif(ns,min = min.value,max = max.value))
    yita0 <- cbind(1,x0)%*%betav
    
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
      prob0 <- 1-pgev(-yita0,loc = 0,scale = 1,shape = xi)
    }
    
    if(model == 'splogit')
    {
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
    
    init <- init.args$init
    logit.fit <- glm(y0~x0,family = binomial(link='logit'))
    probit.fit <- glm(y0 ~ x0, family=binomial(link='probit'))
    robit.fit <- robit.pxem(y0 = y0,x0 = x0,beta0 = init,nu0 = init.args$nu0,tol = 1e-3) 
    
    gev.fit <- gev.mle.xi(x = x0,y = y0,par0 = c(init,init.args$xi0))
    gev.fit.new<- gev.mle.new(y0 = y0,x0 = x0,par0 = c(init,init.args$xi0),maxeval = 3000)
    
    splogit.fit <- splogit.mle(y0 = y0,x0 = x0,par0 = init,intervalr = init.args$intervalr)
   
    bs.nc <- nknots+deg-1
    delta0 <- rep(1/bs.nc,bs.nc)
    
    pspline.lam<- spline.aic(y0,x0,deg=deg,nknots = nknots,kp=1e6,lam.interval = c(0,10000) ,delta0 = delta0,boundary = range(x0),toll =1e-4,  monotone = TRUE)
    
    pspline.fit <- splinelink(y0,x0,deg=deg,nknots = nknots,kp=1e6,lambda = pspline.lam ,delta0 = delta0,boundary = range(x0), monotone = TRUE)
    
   # bgeva.fun <- function(tauv)
   # {
   #   bgeva.fit <- bgeva(y0 ~ s(x0),tau = tauv)
   #   return(bgeva.fit$logL)
   # }

    fb <- optimize(f = bgeva.fun,interval = c(-2,1))
    bgeva.fit <- bgeva(y0~ s(x0),tau = fb$minimum)
    
    prob.gev <-  pgev(q = bgeva.fit$eta,loc = 0,scale = 1,shape = fb$minimum)
    
  
    boundary1.n1[s] <-  min(1-gev.fit$est[3]*cbind(1,x0)%*%gev.fit$est[1:2])
    boundary2.n1[s] <- min(1-gev.fit.new$est[3]*cbind(1,x0)%*%gev.fit.new$est[1:2])
    
    mse.logit[s] <- mean((logit.fit$fitted.values - prob0)^2)
    mse.probit[s] <- mean((probit.fit$fitted.values - prob0)^2)
    mse.robit[s] <- mean((robit.fit$fitted.values - prob0)^2)
    mse.gev[s] <- mean((gev.fit$fitted.values - prob0)^2)
    mse.gev.new[s] <- mean((gev.fit.new$fitted.values - prob0)^2)
    mse.splogit[s] <- mean((splogit.fit$fitted.values - prob0)^2)
    mse.pspline[s] <- mean((pspline.fit$fitted.values - prob0)^2)
   
   aic.logit[s] <- logit.fit$aic
   aic.probit[s] <- probit.fit$aic
   aic.robit[s] <- robit.fit$aic
   aic.gev[s] <- gev.fit$aic
   aic.gev.new[s] <- gev.fit.new$aic
   aic.splogit[s] <- splogit.fit$aic
   aic.pspline[s] <- AIC.value( pspline.fit,y0)
   
    
    
    max.logit[s] <- max(abs(logit.fit$fitted.values - prob0))
    max.probit[s] <- max(abs(probit.fit$fitted.values - prob0))
    max.robit[s] <- max(abs(robit.fit$fitted.values - prob0))
    max.gev[s] <- max(abs(gev.fit$fitted.values - prob0))
    max.gev.new[s] <- max(abs(gev.fit.new$fitted.values - prob0))
    max.splogit[s] <- max(abs(splogit.fit$fitted.values - prob0))
    max.pspline[s] <- max(abs(pspline.fit$fitted.values - prob0))
    
    ks.logit[s] <- ks.test(prob0,logit.fit$fitted.values)$p.value
    ks.probit[s] <- ks.test(prob0,probit.fit$fitted.values)$p.value
    ks.robit[s] <- ks.test(prob0,robit.fit$fitted.values)$p.value
    ks.gev[s] <- ks.test(prob0,gev.fit$fitted.values)$p.value
    ks.gev.new[s] <- ks.test(prob0,gev.fit.new$fitted.values)$p.value
    ks.splogit[s] <- ks.test(prob0,splogit.fit$fitted.values)$p.value
    ks.pspline[s] <- ks.test(prob0,pspline.fit$fitted.values)$p.value
    
    xmat <- cbind(1,x0)
    qtle.logit[s,] <- quantile(xmat%*%coef(logit.fit),probs = quantile)
    qtle.probit[s,] <- quantile(xmat%*%coef(probit.fit),probs = quantile)
    qtle.robit[s,] <- quantile(robit.fit$eta,probs = quantile)
    qtle.gev[s,] <- quantile(gev.fit$eta,probs = quantile)
    qtle.gev.new[s,] <- quantile(gev.fit.new$eta,probs = quantile)
    qtle.splogit[s,] <- quantile(splogit.fit$eta,probs = quantile)
    qtle.pspline[s,] <- quantile(pspline.fit$eta,probs = quantile)
    
    
    grv1.n1[s,] <- gev.fit$gr
    grv2.n1[s,] <- gev.fit.new$gr
    splogit.rv.n1[s,] <- splogit.fit$gr
    
    print(s)
  }
  
  mse.mat <- cbind(mse.logit,mse.probit,mse.robit,
                   mse.gev,mse.gev.new,mse.splogit,mse.pspline)
  max.mat <- cbind(max.logit,max.probit,max.robit,
                   max.gev,max.gev.new,max.splogit,max.pspline)
  
  pmat <- cbind(ks.logit,ks.probit,ks.robit,
                ks.gev,ks.gev.new,ks.splogit,ks.pspline)
  
  qtle.list <- list(logit = qtle.logit,probit = qtle.probit,robit = qtle.robit,gev = qtle.gev,gev.new = qtle.gev.new,splogit = qtle.splogit,pspline = qtle.pspline)
  outls <- list(mse.mat = mse.mat, max.mat = max.mat,pmat = pmat, 
                gr1 = grv1.n1, gr2 = grv2.n1, boundary1 = boundary1.n1, boundary2 = boundary2.n1,splogit.rv = splogit.rv.n1,qtle=qtle.list)
  return(outls)
}


