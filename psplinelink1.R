#### algorithm #### 
#### spline not pspline, 5 knots ###
library(actuar)
library(splines)
library(MASS)
## qv is value of quantile value of qv ####
psplinelink1<- function(y0,xmat,qv=1,deg = 3,d.value =5,kp = 1e6,nknots=5,monotone=TRUE,beta0,delta0,tol = 1e-8,lambda=2000)
{
  xmats <- scale(xmat) 
  atu <- quantile(sqrt(rowSums(xmats^2)),qv)
 
  bs.nc <- nknots+deg-1
  delta0 <- rep(0,bs.nc)
  knots <- seq(0,1,length.out = nknots)
  
  beta.old <- beta0/sqrt(sum(beta0^2))
  delta.old <- delta0
  
  ##### link function #### 
  repeat
  {   
    beta.comp <- beta.old
    delta.comp <- delta.old
    
    eta.old <- xmats%*%beta.comp
    q.old <- (eta.old/atu+1)/2
    Ut <- pgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )
    
    bs0 <- bs(Ut,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots = c(0,1),intercept=TRUE)
    
    if(!monotone)
    {
      repeat
      {       
        bs.eta <- bs0%*%delta.old
        bs.mu <- exp(bs.eta)/(1+exp(bs.eta))
        index <-  bs.mu > .Machine$double.eps
        wt <- diag(as.numeric(bs.mu[index]*(1-bs.mu[index])))
        z <- bs.eta[index]+(y0[index]-bs.mu[index])/as.numeric(bs.mu[index]*(1-bs.mu[index])) 
        delta.update <- ginv(t(bs0[index,])%*%wt%*%bs0[index,])%*%t(bs0[index,])%*%wt%*%z
        diff.value <- sqrt(sum((delta.update-delta.old)^2))
        delta.old <- delta.update
        if(diff.value <= tol){break}
      }  
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
      
      repeat
      {
        bs.eta <- bs0%*%delta.old
        bs.mu <- exp(bs.eta)/(1+exp(bs.eta))
        wt <- diag(as.numeric(bs.mu*(1-bs.mu)))
        z <- bs.eta+(y0-bs.mu)/as.numeric(bs.mu*(1-bs.mu)) 
        Vmat <- diag(as.numeric(Dmat%*%delta.old) <0 ) 
        delta.update <- ginv(t(bs0)%*%wt%*%bs0 + lambda*t(Dmat1)%*%(Dmat1) +kp*t(Dmat)%*%Vmat%*%Dmat)%*%t(bs0)%*%wt%*%z
        diff.value <- sqrt(sum((delta.update-delta.old)^2))
        delta.old <- delta.update
        if(diff.value <= tol){break}
        
      }
      
    }
    
    ### update beta ####
    
    repeat
    {
      eta.old<- xmats%*%beta.old
      eta.stand <- eta.old
      q.old <- (eta.stand/atu+1)/2
      Ut <- pgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )   
      bs.old <- bs(Ut,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots =c(0,1),intercept=TRUE)
      
      muhat <- exp(bs.old %*% delta.update)/(1+exp(bs.old %*% delta.update))
      bs.deriv <- splineDesign(knots= c(rep(0,4),knots[c(-1,-length(knots))],rep(1,4)), Ut, ord = 4, derivs=rep(1,length(y0)),outer.ok=TRUE)
      fun.deriv <- bs.deriv%*% delta.update
      den.value<- dgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )/(2*atu)
      du.eta <- muhat*(1-muhat)* fun.deriv*den.value
      z.beta <- eta.stand +(y0 - muhat)/du.eta
      wt.beta <- diag(as.numeric(muhat*(1-muhat)*fun.deriv^2*den.value^2))
      beta.update <- solve(t(xmats)%*%wt.beta%*%xmats)%*%t(xmats)%*%wt.beta%*%z.beta
      beta.update <- beta.update/sqrt(sum(beta.update^2))
      diff.value <- mean((beta.update-beta.old)^2)
      beta.old <- beta.update
      if(diff.value <= tol){break} 
    }
    
    diff.total<- sqrt(sum((beta.update - beta.comp)^2) + sum((delta.update - delta.comp)^2))
    if(diff.total<= tol){break} 
  }
  
  eta <- xmats%*%beta.update
  fitted.values <- muhat
  res <- list(eta= eta,est = beta.update,fitted.values = fitted.values)
  return(res)
}

#### pspline using GCV ###
pspline.gcv <- function(y0,xmat,qv=1,deg = 3,d.value =5,nknots=5,kp=1e6,monotone=TRUE,beta0,delta0,tol = 1e-8,lam.interval)
{  
  bs.nc <- nknots+deg-1
  xmats <- scale(xmat)
  atu <- quantile(sqrt(rowSums(xmats^2)),qv)
  delta0 <- rep(0,bs.nc)
  knots <- seq(0,1,length.out = nknots)
  
  lam.fun<- function(lambda)
  {   
    beta.old <- beta0/sqrt(sum(beta0^2))
    delta.old <- delta0
    
    ##### link function #### 
    repeat
    {   
      beta.comp <- beta.old
      delta.comp <- delta.old
      
      eta.old <- xmats%*%beta.comp
      q.old <- (eta.old/atu+1)/2
      Ut <- pgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )
      
      bs0 <- bs(Ut,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots = c(0,1),intercept=TRUE)
      
      if(!monotone)
      {
        repeat
        {       
          bs.eta <- bs0%*%delta.old
          bs.mu <- exp(bs.eta)/(1+exp(bs.eta))
          index <-  bs.mu > .Machine$double.eps
          wt <- diag(as.numeric(bs.mu[index]*(1-bs.mu[index])))
          z <- bs.eta[index]+(y0[index]-bs.mu[index])/as.numeric(bs.mu[index]*(1-bs.mu[index])) 
          delta.update <- ginv(t(bs0[index,])%*%wt%*%bs0[index,])%*%t(bs0[index,])%*%wt%*%z
          diff.value <- sqrt(sum((delta.update-delta.old)^2))
          delta.old <- delta.update
          if(diff.value <= tol){break}
        }  
        
        Hmat <- solve(t(bs0)%*%wt%*%bs0 + lambda*t(Dmat1)%*%(Dmat1))%*%t(bs0)%*%wt%*%bs0
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
        
        repeat
        {
          bs.eta <- bs0%*%delta.old
          bs.mu <- exp(bs.eta)/(1+exp(bs.eta))
          wt <- diag(as.numeric(bs.mu*(1-bs.mu)))
          z <- bs.eta+(y0-bs.mu)/as.numeric(bs.mu*(1-bs.mu)) 
          Vmat <- diag(as.numeric(Dmat%*%delta.old) <0 ) 
          delta.update <- ginv(t(bs0)%*%wt%*%bs0 + lambda*t(Dmat1)%*%(Dmat1) +kp*t(Dmat)%*%Vmat%*%Dmat)%*%t(bs0)%*%wt%*%z
          diff.value <- sqrt(sum((delta.update-delta.old)^2))
          delta.old <- delta.update
          if(diff.value <= tol){break}
          
        }
        
        Hmat <- solve(t(bs0)%*%wt%*%bs0 + lambda*t(Dmat1)%*%(Dmat1)+kp*t(Dmat)%*%Vmat%*%Dmat)%*%t(bs0)%*%wt%*%bs0
        traceH <- sum(diag(Hmat))
      }
      
      ### update beta ####
      
      repeat
      {
        eta.old<- xmats%*%beta.old
        eta.stand <- eta.old
        q.old <- (eta.stand/atu+1)/2
        Ut <- pgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )   
        bs.old <- bs(Ut,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots =c(0,1),intercept=TRUE)
        
        muhat <- exp(bs.old %*% delta.update)/(1+exp(bs.old %*% delta.update))
        bs.deriv <- splineDesign(knots= c(rep(0,4),knots[c(-1,-length(knots))],rep(1,4)), Ut, ord = 4, derivs=rep(1,length(y0)),outer.ok=TRUE)
        fun.deriv <- bs.deriv%*% delta.update
        den.value<- dgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )/(2*atu)
        du.eta <- muhat*(1-muhat)* fun.deriv*den.value
        z.beta <- eta.stand +(y0 - muhat)/du.eta
        wt.beta <- diag(as.numeric(muhat*(1-muhat)*fun.deriv^2*den.value^2))
        beta.update <- solve(t(xmats)%*%wt.beta%*%xmats)%*%t(xmats)%*%wt.beta%*%z.beta
        beta.update <- beta.update/sqrt(sum(beta.update^2))
        diff.value <- mean((beta.update-beta.old)^2)
        beta.old <- beta.update
        if(diff.value <= tol){break} 
      }
      
      diff.total<- sqrt(sum((beta.update - beta.comp)^2) + sum((delta.update - delta.comp)^2))
      if(diff.total<= tol){break} 
    }
    
    eta <- xmats%*%beta.update
    fitted.values <- muhat
    gcv <- mean((y0 - muhat)^2)/(1-traceH/length(y0))^2
    return(gcv)
  }
  
  out.gcv <- optimize(f = lam.fun,interval = lam.interval)
  lam.value <- out.gcv$minimum
  return(lam.value)
}

print(c('psplinelink1','pspline.gcv'))



