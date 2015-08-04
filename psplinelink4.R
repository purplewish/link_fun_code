###############  the algorithm in the paper ######

#### algorithm #### 
library(actuar)
library(splines)
library(MASS)
## qv is value of quantile value of qv ####
# cat is the name of categorical data ###
psplinelink4<- function(y0,xmat0,deg = 3,kp = 1e6,nknots=10,
                        monotone=TRUE,beta0,delta0,
                        tol = 1e-8,lambda=20,MaxIter=1000,boundary)
{
  
  bs.nc <- nknots+deg-1
  delta0 <- rep(0,bs.nc)

  beta.old <- beta0
  delta.old <- delta0
  
  xmat <- cbind(1,as.matrix(xmat0))
  nr <- length(y0)
  
  ##### link function #### 
  for(j in 1:MaxIter)
  {
    eta.old <- xmat%*%beta.old
    eta.std <- as.numeric(scale(eta.old)*sqrt(nr/(nr-1)))
    knots <- seq(min(eta.std),max(eta.std),length.out = nknots)
    bs.old <- bs(eta.std,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots = range(eta.std) ,intercept=TRUE)
    
    for(j1 in 1:MaxIter)
    {
      if(!monotone)
      {
        Dmat1 <- matrix(0,ncol = bs.nc, nrow = bs.nc -2)
        Dmat1[cbind(1:(bs.nc-2),1:(bs.nc-2))] <- 1
        Dmat1[cbind(1:(bs.nc-2),2:(bs.nc-1))] <- -2
        Dmat1[cbind(1:(bs.nc-2),3:bs.nc)] <- 1
        
        bs.eta <- bs.old%*%delta.old
        bs.mu <- exp(bs.eta)/(1+exp(bs.eta))
        index <-  bs.mu > .Machine$double.eps
        wt <- diag(as.numeric(bs.mu[index]*(1-bs.mu[index])))
        z <- bs.eta[index]+(y0[index]-bs.mu[index])/as.numeric(bs.mu[index]*(1-bs.mu[index])) 
        delta.update <- ginv(t(bs.old[index,])%*%wt%*%bs.old[index,]+lambda*t(Dmat1)%*%(Dmat1))%*%t(bs.old[index,])%*%wt%*%z
        Hmat <- solve(t(bs.old)%*%wt%*%bs.old + lambda*t(Dmat1)%*%(Dmat1))%*%t(bs.old)%*%wt%*%bs.old
        traceH <- sum(diag(Hmat))
      }
      
      if(monotone)
      {
        Dmat <- matrix(0,ncol=bs.nc,nrow=bs.nc-1)
        diag(Dmat) <- -1
        Dmat[cbind(1:(bs.nc-1),2:bs.nc)] <- 1
        
        Dmat1 <- matrix(0,ncol = bs.nc, nrow = bs.nc -2)
        Dmat1[cbind(1:(bs.nc-2),1:(bs.nc-2))] <- 1
        Dmat1[cbind(1:(bs.nc-2),2:(bs.nc-1))] <- -2
        Dmat1[cbind(1:(bs.nc-2),3:bs.nc)] <- 1
        
        
        bs.eta <- bs.old%*%delta.old
        bs.mu <- exp(bs.eta)/(1+exp(bs.eta))
        wt <- diag(as.numeric(bs.mu*(1-bs.mu)))
        z <- bs.eta+(y0-bs.mu)/as.numeric(bs.mu*(1-bs.mu)) 
        Vmat <- diag(as.numeric(Dmat%*%delta.old) <0 ) 
        #         cz <- chol(t(bs.old)%*%wt%*%bs.old + lambda*t(Dmat1)%*%(Dmat1) +kp*t(Dmat)%*%Vmat%*%Dmat)
        #         delta.update <- chol2inv(cz)%*%t(bs.old)%*%wt%*%z
        delta.update <- ginv(t(bs.old)%*%wt%*%bs.old + lambda*t(Dmat1)%*%(Dmat1) +kp*t(Dmat)%*%Vmat%*%Dmat)%*%t(bs.old)%*%wt%*%z
        
        Hmat <- solve(t(bs.old)%*%wt%*%bs.old + lambda*t(Dmat1)%*%(Dmat1)+
                        kp*t(Dmat)%*%Vmat%*%Dmat)%*%t(bs.old)%*%wt%*%bs.old
        traceH <- sum(diag(Hmat))
        
      }
      
      diff.delta <- sqrt(sum((delta.update-delta.old)^2))
      if(diff.delta <- tol){break}
      delta.old <- delta.update
      if(j == MaxIter){print('MaxIter reached without convergence')}
    }

    
    ### update beta ####  
    
    muhat <- exp(bs.old %*% delta.update)/(1+exp(bs.old %*% delta.update))
    bs.deriv <- splineDesign(knots= c(rep(min(eta.std),4),knots[c(-1,-length(knots))],rep(max(eta.std),4)), eta.std , ord = 4, derivs=rep(1,nr),outer.ok=TRUE)
    fun.deriv <- bs.deriv%*% delta.update
    du.eta <- muhat*(1-muhat)* fun.deriv
    z.beta <- eta.std +(y0 - muhat)/du.eta
    wt.beta <- diag(as.numeric(muhat*(1-muhat)*fun.deriv^2))
    beta.update <- solve(t(xmat)%*%wt.beta%*%xmat)%*%t(xmat)%*%wt.beta%*%z.beta
    
    
    diff.total<- sqrt(sum((beta.update - beta.old)^2) + sum((delta.update - delta.old)^2))
    if(diff.total<= tol){break} 
    beta.old <- beta.update
    delta.old <- delta.update
    if(j == MaxIter){print('MaxIter reached without convergence')}
    
  }
  
  
  eta <- xmat%*%beta.update
  eta.std <- as.numeric(scale(eta)*sqrt(nr/(nr-1)))
  eta.center <- mean(eta)
  eta.sds <- sqrt((nr-1)/nr)*sd(eta)
  bs0 <- bs(eta.std,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots = range(eta.std),intercept=TRUE)
  fitted.values <- exp(bs0%*%delta.update)/(1+exp(bs0%*%delta.update))
  
  res <- list(eta= eta,est = beta.update,delta=delta.update,fitted.values = fitted.values,
              deg=deg,kp=kp,nknots=nknots,knots=knots,traceH=traceH,eta.center=eta.center,eta.sds=eta.sds)
  
  return(res)
}

#### pspline using GCV ###
pspline.gcv4 <- function(y0,xmat0,deg = 3,nknots=5,kp=1e6,
 monotone=TRUE,beta0,delta0,tol = 1e-8,lamv=seq(5,100,length.out = 20),MaxIter = 100)
{  
  
  bs.nc <- nknots+deg-1
  delta0 <- rep(0,bs.nc)
  
  xmat <- cbind(1,as.matrix(xmat0))
  nr <- length(y0)
  
  lam.fun<- function(lambda)
  {   
    beta.old <- beta0
    delta.old <- delta0
  
    ##### link function #### 

    for(j in 1:MaxIter)
    {
      eta.old <- xmat%*%beta.old
      eta.std <- as.numeric(scale(eta.old)*sqrt(nr/(nr-1)))
      knots <- seq(min(eta.std),max(eta.std),length.out = nknots)
      bs.old <- bs(eta.std,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots = range(eta.std) ,intercept=TRUE)
      
      for(j1 in 1:MaxIter)
      {
        if(!monotone)
        {
          Dmat1 <- matrix(0,ncol = bs.nc, nrow = bs.nc -2)
          Dmat1[cbind(1:(bs.nc-2),1:(bs.nc-2))] <- 1
          Dmat1[cbind(1:(bs.nc-2),2:(bs.nc-1))] <- -2
          Dmat1[cbind(1:(bs.nc-2),3:bs.nc)] <- 1
          
          bs.eta <- bs.old%*%delta.old
          bs.mu <- exp(bs.eta)/(1+exp(bs.eta))
          index <-  bs.mu > .Machine$double.eps
          wt <- diag(as.numeric(bs.mu[index]*(1-bs.mu[index])))
          z <- bs.eta[index]+(y0[index]-bs.mu[index])/as.numeric(bs.mu[index]*(1-bs.mu[index])) 
          delta.update <- ginv(t(bs.old[index,])%*%wt%*%bs.old[index,]+lambda*t(Dmat1)%*%(Dmat1))%*%t(bs.old[index,])%*%wt%*%z
          Hmat <- solve(t(bs.old)%*%wt%*%bs.old + lambda*t(Dmat1)%*%(Dmat1))%*%t(bs.old)%*%wt%*%bs.old
          traceH <- sum(diag(Hmat))
        }
        
        if(monotone)
        {
          Dmat <- matrix(0,ncol=bs.nc,nrow=bs.nc-1)
          diag(Dmat) <- -1
          Dmat[cbind(1:(bs.nc-1),2:bs.nc)] <- 1
          
          Dmat1 <- matrix(0,ncol = bs.nc, nrow = bs.nc -2)
          Dmat1[cbind(1:(bs.nc-2),1:(bs.nc-2))] <- 1
          Dmat1[cbind(1:(bs.nc-2),2:(bs.nc-1))] <- -2
          Dmat1[cbind(1:(bs.nc-2),3:bs.nc)] <- 1
          
          
          bs.eta <- bs.old%*%delta.old
          bs.mu <- exp(bs.eta)/(1+exp(bs.eta))
          wt <- diag(as.numeric(bs.mu*(1-bs.mu)))
          z <- bs.eta+(y0-bs.mu)/as.numeric(bs.mu*(1-bs.mu)) 
          Vmat <- diag(as.numeric(Dmat%*%delta.old) <0 ) 
          #         cz <- chol(t(bs.old)%*%wt%*%bs.old + lambda*t(Dmat1)%*%(Dmat1) +kp*t(Dmat)%*%Vmat%*%Dmat)
          #         delta.update <- chol2inv(cz)%*%t(bs.old)%*%wt%*%z
          delta.update <- ginv(t(bs.old)%*%wt%*%bs.old + lambda*t(Dmat1)%*%(Dmat1) +kp*t(Dmat)%*%Vmat%*%Dmat)%*%t(bs.old)%*%wt%*%z
          
          Hmat <- solve(t(bs.old)%*%wt%*%bs.old + lambda*t(Dmat1)%*%(Dmat1)+
                          kp*t(Dmat)%*%Vmat%*%Dmat)%*%t(bs.old)%*%wt%*%bs.old
          traceH <- sum(diag(Hmat))
          
        }
        
        diff.delta <- sqrt(sum((delta.update-delta.old)^2))
        if(diff.delta <- tol){break}
        delta.old <- delta.update
        if(j == MaxIter){print('MaxIter reached without convergence')}
      }
      
      
      ### update beta ####  
      
      muhat <- exp(bs.old %*% delta.update)/(1+exp(bs.old %*% delta.update))
      bs.deriv <- splineDesign(knots= c(rep(min(eta.std),4),knots[c(-1,-length(knots))],rep(max(eta.std),4)), eta.std , ord = 4, derivs=rep(1,nr),outer.ok=TRUE)
      fun.deriv <- bs.deriv%*% delta.update
      du.eta <- muhat*(1-muhat)* fun.deriv
      z.beta <- eta.std +(y0 - muhat)/du.eta
      wt.beta <- diag(as.numeric(muhat*(1-muhat)*fun.deriv^2))
      beta.update <- solve(t(xmat)%*%wt.beta%*%xmat)%*%t(xmat)%*%wt.beta%*%z.beta
      
      
      diff.total<- sqrt(sum((beta.update - beta.old)^2) + sum((delta.update - delta.old)^2))
      if(diff.total<= tol){break} 
      beta.old <- beta.update
      delta.old <- delta.update
      if(j == MaxIter){print('MaxIter reached without convergence')}
      
    }
    eta <- xmat%*%beta.update
    eta.std <- as.numeric(scale(eta)*sqrt(nr/(nr-1)))
    bs0 <- bs(eta.std,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots = range(eta.std),intercept=TRUE)
    fitted.values <- exp(bs0%*%delta.update)/(1+exp(bs0%*%delta.update))
    gcv <- mean((y0 - fitted.values)^2)/(1-traceH/nr)^2
    return(gcv)
  }
  
  
  gcvv <- unlist(lapply(lamv,lam.fun))
  #out.gcv <- optimize(f = lam.fun,interval = lam.interval)
  #lam.value <- out.gcv$minimum
  lam.value <- lamv[which.min(gcvv)]
  
  return(lam.value)
}


##### predict based on pspline link function ######

predict.pspline4 <- function(est.obj,newdata)
{
  deg <- est.obj$deg
  nknots <- est.obj$nknots
  knots <- est.obj$knots
  est <- est.obj$est
  delta.est <- est.obj$delta
  
  eta.center <- est.obj$eta.center
  eta.sds <- est.obj$eta.sds
  
  eta <- cbind(1,newdata)%*%est
  eta.std <- (eta - eta.center)/eta.sds
  bs.value <- bs(eta.std,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots = knots[c(1,length(knots))],intercept=TRUE)
  eta.bs <- bs.value%*%delta.est
  prob.est <- exp(eta.bs)/(1+exp(eta.bs))
  return(prob.est)
  
}

print(c('psplinelink4','pspline.gcv4','predict.pspline4'))