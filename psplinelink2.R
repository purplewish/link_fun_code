#### algorithm #### 
library(actuar)
library(splines)
library(MASS)
## qv is value of quantile value of qv ####
# cat is the name of categorical data ###
psplinelink2<- function(y0,xmat,qv=1,deg = 3,kp = 1e6,nknots=5,
                        monotone=TRUE,beta0,delta0,catv=NULL,
                        tol = 1e-8,lambda=20,MaxIter=1000)
{
  xmats <- scale(xmat) 
  d.value <- ncol(xmat)
  xmats[,catv] <- xmat[,catv]
  center <- attr(xmats,'scaled:center')
  sd.value <- attr(xmats,'scaled:scale')
  atu <- quantile(sqrt(rowSums(xmats^2)),0.95)
  
  bs.nc <- nknots+deg-1
  delta0 <- rep(0,bs.nc)
  knots <- seq(0,1,length.out = nknots)
  
  beta.old <- beta0/sqrt(sum(beta0^2))
  delta.old <- delta0
  
  ##### link function #### 
  for(j in 1:MaxIter)
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
      
      for(j1 in 1:MaxIter)
      {
        
        bs.eta <- bs0%*%delta.old
        bs.mu <- exp(bs.eta)/(1+exp(bs.eta))
        wt <- diag(as.numeric(bs.mu*(1-bs.mu)))
        z <- bs.eta+(y0-bs.mu)/as.numeric(bs.mu*(1-bs.mu)) 
        Vmat <- diag(as.numeric(Dmat%*%delta.old) <0 ) 
#         cz <- chol(t(bs0)%*%wt%*%bs0 + lambda*t(Dmat1)%*%(Dmat1) +kp*t(Dmat)%*%Vmat%*%Dmat)
#         delta.update <- chol2inv(cz)%*%t(bs0)%*%wt%*%z
     delta.update <- ginv(t(bs0)%*%wt%*%bs0 + lambda*t(Dmat1)%*%(Dmat1) +kp*t(Dmat)%*%Vmat%*%Dmat)%*%t(bs0)%*%wt%*%z
        
        
        diff.value <- sqrt(sum((delta.update-delta.old)^2))
        delta.old <- delta.update
print(delta.update)
        if(diff.value <= tol){break}
        if(j1 == MaxIter){print('MaxIter reached without convergence')}
      }
      
    }
    
    ### update beta ####

#    gamma.fun <- function(gamma)
#    {
#      for(j2 in 1:MaxIter)
#      {
#        eta.old<- xmats%*%beta.old
#        eta.stand <- eta.old
#        q.old <- (eta.stand/atu+1)/2
#        Ut <- pgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )   
#        bs.old <- bs(Ut,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots =c(0,1),intercept=TRUE)
#        
#        muhat <- exp(bs.old %*% delta.update)/(1+exp(bs.old %*% delta.update))
#        bs.deriv <- splineDesign(knots= c(rep(0,4),knots[c(-1,-length(knots))],rep(1,4)), Ut, ord = 4, derivs=rep(1,length(y0)),outer.ok=TRUE)
#        fun.deriv <- bs.deriv%*% delta.update
#        den.value<- dgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )/(2*atu)
#        du.eta <- muhat*(1-muhat)* fun.deriv*den.value
#        z.beta <- eta.stand +(y0 - muhat)/du.eta
#        wt.beta <- diag(as.numeric(muhat*(1-muhat)*fun.deriv^2*den.value^2))
#        # beta.update <- solve(t(xmats)%*%wt.beta%*%xmats)%*%t(xmats)%*%wt.beta%*%z.beta
#        beta.update<- solve(t(xmats)%*%wt.beta%*%xmats - 2*gamma*diag(1,2))%*%t(xmats)%*%wt.beta%*%z.beta
#        beta.update <- beta.update/sqrt(sum(beta.update^2))
#        diff.value <- sqrt(sum((beta.update-beta.old)^2))
#        beta.old <- beta.update
#      }
#  
#    }


   beta.fun <- function(beta.update)
   {
    
     eta.old<- xmats%*%beta.update
     eta.stand <- eta.old
     q.old <- (eta.stand/atu+1)/2
     Ut <- pgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )   
     bs.old <- bs(Ut,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots =c(0,1),intercept=TRUE)
     
     muhat <- exp(bs.old %*% delta.update)/(1+exp(bs.old %*% delta.update))
     
     nllv <- sum(-y0*log(muhat) - (1-y0)*log(1-muhat))
     
     bs.deriv <- splineDesign(knots= c(rep(0,4),knots[c(-1,-length(knots))],rep(1,4)), Ut, ord = 4, derivs=rep(1,length(y0)),outer.ok=TRUE)
     fun.deriv <- bs.deriv%*% delta.update
     den.value<- dgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )/(2*atu)
     du.eta <- muhat*(1-muhat)* fun.deriv*den.value
     z.beta <- (y0 - muhat)/du.eta
     wt.beta <- diag(as.numeric(muhat*(1-muhat)*fun.deriv^2*den.value^2))
     gr.value <- -t(xmats)%*%wt.beta%*%z.beta
     
     return( list("objective"=nllv, 
                  "gradient"= gr.value ) )

   }



    eval_g1 <- function(x) 
    {
    return( list("constraints"= sum(x^2) - 1,
               "jacobian"= 2*x) )
     }

  res <- nloptr( x0=beta.old, eval_f=beta.fun, eval_g_eq = eval_g1, 
               opts = list("algorithm"="NLOPT_LD_SLSQP", "check_derivatives"=TRUE,maxeval=maxeval)) 

 
  beta.update <- res$solution
  beta.old <- beta.update
 
    diff.total<- sqrt(sum((beta.update - beta.comp)^2) + sum((delta.update - delta.comp)^2))
    if(diff.total<= tol){break} 
    if(j == MaxIter){print('MaxIter reached without convergence')}    
  }
  
  
  eta <- xmats%*%beta.update
  q.value <- (eta/atu+1)/2
  Ut <- pgenbeta(q.value,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )
  bs0 <- bs(Ut,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots = c(0,1),intercept=TRUE)
  fitted.values <- exp(bs0%*%delta.update)/(1+exp(bs0%*%delta.update))
  res <- list(eta= eta,est = beta.update,delta=delta.update,fitted.values = fitted.values,
              qv=qv,deg=deg,kp=kp,nknots=nknots,knots=knots,d.value=d.value,atu=atu,
              center=center,sd.value = sd.value,catv=catv)
  
  return(res)
}

#### pspline using GCV ###
pspline.gcv2 <- function(y0,xmat,qv=1,deg = 3,nknots=5,kp=1e6,catv=NULL,
                        monotone=TRUE,beta0,delta0,tol = 1e-8,lamv=seq(5,100,length.out = 20),MaxIter = 100)
{  
  xmats <- scale(xmat) 
  xmats[,catv] <- xmat[,catv]
  d.value <- ncol(xmat)
  atu <- quantile(sqrt(rowSums(xmats^2)),qv)
  bs.nc <- nknots+deg-1
  delta0 <- rep(0,bs.nc)
  knots <- seq(0,1,length.out = nknots)
  
  lam.fun<- function(lambda)
  {   
    beta.old <- beta0/sqrt(sum(beta0^2))
    delta.old <- delta0
    
    ##### link function #### 
    for(j in 1:MaxIter)
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
         traceH<- sum(diag(Hmat))
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
        
        for(j1 in 1:MaxIter)
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
          if(j1 == MaxIter){print('MaxIter reached without convergence')}
        }
        
        
        Hmat <- solve(t(bs0)%*%wt%*%bs0 + lambda*t(Dmat1)%*%(Dmat1)+kp*t(Dmat)%*%Vmat%*%Dmat)%*%t(bs0)%*%wt%*%bs0
        traceH <- sum(diag(Hmat))
      }
      
      ### update beta ####

      beta.fun <- function(beta.update)
      {
        
        eta.old<- xmats%*%beta.update
        eta.stand <- eta.old
        q.old <- (eta.stand/atu+1)/2
        Ut <- pgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )   
        bs.old <- bs(Ut,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots =c(0,1),intercept=TRUE)
        
        muhat <- exp(bs.old %*% delta.update)/(1+exp(bs.old %*% delta.update))
        
        nllv <- sum(-y0*log(muhat) - (1-y0)*log(1-muhat))
        
        bs.deriv <- splineDesign(knots= c(rep(0,4),knots[c(-1,-length(knots))],rep(1,4)), Ut, ord = 4, derivs=rep(1,length(y0)),outer.ok=TRUE)
        fun.deriv <- bs.deriv%*% delta.update
        den.value<- dgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )/(2*atu)
        du.eta <- muhat*(1-muhat)* fun.deriv*den.value
        z.beta <- (y0 - muhat)/du.eta
        wt.beta <- diag(as.numeric(muhat*(1-muhat)*fun.deriv^2*den.value^2))
        gr.value <- -t(xmats)%*%wt.beta%*%z.beta
        
        return( list("objective"=nllv, 
                     "gradient"= gr.value ) )
        
      }
      
      
      
      eval_g1 <- function(x) 
      {
        return( list("constraints"= sum(x^2) - 1,
                     "jacobian"= 2*x) )
      }
      
      res <- nloptr( x0=beta.old, eval_f=beta.fun, eval_g_eq = eval_g1, 
                     opts = list("algorithm"="NLOPT_LD_SLSQP", "check_derivatives"=TRUE,maxeval=maxeval)) 
      
      
      beta.update <- res$solution
      beta.old <- beta.update
      
      
      diff.total<- sqrt(sum((beta.update - beta.comp)^2) + sum((delta.update - delta.comp)^2))
      if(diff.total<= tol){break} 
      if(j == MaxIter){print('MaxIter reached without convergence')}
      
    }
    
    
    eta <- xmats%*%beta.update
    fitted.values <- muhat
    gcv <- mean((y0 - muhat)^2)/(1-traceH/length(y0))^2
    return(gcv)
  }
  
  
  gcvv <- unlist(lapply(lamv,lam.fun))
  #out.gcv <- optimize(f = lam.fun,interval = lam.interval)
  #lam.value <- out.gcv$minimum
  lam.value <- lamv[which.min(gcvv)]
  
  return(lam.value)
}




##### predict based on pspline link function ######

predict.pspline1 <- function(est.obj,newdata)
{
  
  deg <- est.obj$deg
  nknots <- est.obj$nknots
  knots <- est.obj$knots
  d.value <- est.obj$d.value
  qv <- est.obj$qv 
  est <- est.obj$est
  atu <- est.obj$atu
  delta.est <- est.obj$delta
  catv<- est.obj$catv
  center <- est.obj$center
  sd.value <- est.obj$sd.value
  
  newdatas <- t((t(newdata) - center)/sd.value)
  newdatas[,catv] <- newdata[,catv]
  
  eta.new <- newdatas%*%est
  q.value <- (eta.new/atu+1)/2
  Ut <- pgenbeta(q.value,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )
  
  bs.value <- bs(Ut,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots = c(0,1),intercept=TRUE)
  
  eta.bs <- bs.value%*%delta.est
  prob.est <- exp(eta.bs)/(1+exp(eta.bs))
  return(prob.est)
  
}

print(c('psplinelink1','pspline.gcv','predict.pspline1'))
