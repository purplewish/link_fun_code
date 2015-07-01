###############  one step update in each step (there are two steps) ######

#### algorithm #### 
library(actuar)
library(splines)
library(MASS)
## qv is value of quantile value of qv ####
# cat is the name of categorical data ###
### size >1 , binomial ###
psplinelink3<- function(y0,xmat,size=1,qv=1,deg = 3,kp = 1e6,nknots=10,
                        monotone=TRUE,beta0,delta0,catv=NULL,
                        tol = 1e-8,lambda=20,MaxIter=1000,dd)
{
  xmats <- scale(xmat) 
  d.value <- ncol(xmat)
  xmats[,catv] <- xmat[,catv]
  center <- attr(xmats,'scaled:center')
  sd.value <- attr(xmats,'scaled:scale')
  atu <- quantile(sqrt(rowSums(xmats^2)),qv)
  
  bs.nc <- nknots+deg-1
  delta0 <- rep(0,bs.nc)
  knots <- seq(0,1,length.out = nknots)
  
  beta.old <- beta0/sqrt(sum(beta0^2))
  delta.old <- delta0
  
  ##### link function #### 
  for(j in 1:MaxIter)
  {
  
    eta.old <- xmats%*%beta.old
    q.old <- (eta.old/atu+1)/2
    Ut <- pgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )
    
    bs.old <- bs(Ut,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots = c(0,1),intercept=TRUE)
    
    if(!monotone)
    {
      Dmat1 <- matrix(0,ncol = bs.nc, nrow = bs.nc -2)
      Dmat1[cbind(1:(bs.nc-2),1:(bs.nc-2))] <- 1
      Dmat1[cbind(1:(bs.nc-2),2:(bs.nc-1))] <- -2
      Dmat1[cbind(1:(bs.nc-2),3:bs.nc)] <- 1
     
        bs.eta <- bs.old%*%delta.old
        bs.mu <- exp(bs.eta)/(1+exp(bs.eta))
        wt <- as.numeric(size*bs.mu*(1-bs.mu))
        z <- bs.eta+(y0-size*bs.mu)/as.numeric(size*bs.mu*(1-bs.mu)) 
        delta.update <- ginv(t(wt*bs.old)%*%bs.old+lambda*t(Dmat1)%*%(Dmat1))%*%t(wt*bs.old)%*%z
        Hmat <- solve(t(wt*bs.old)%*%bs.old + lambda*t(Dmat1)%*%(Dmat1))%*%t(wt*bs.old)%*%bs.old
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
        wt <- as.numeric(size*bs.mu*(1-bs.mu))
        z <- bs.eta+(y0-size*bs.mu)/as.numeric(size*bs.mu*(1-bs.mu)) 
        Vmat <- diag(as.numeric(Dmat%*%delta.old) <0 ) 
#         cz <- chol(t(bs.old)%*%wt%*%bs.old + lambda*t(Dmat1)%*%(Dmat1) +kp*t(Dmat)%*%Vmat%*%Dmat)
#         delta.update <- chol2inv(cz)%*%t(bs.old)%*%wt%*%z
                delta.update <- ginv(t(wt*bs.old)%*%bs.old + lambda*t(Dmat1)%*%(Dmat1) +kp*t(Dmat)%*%Vmat%*%Dmat)%*%t(wt*bs.old)%*%z
        
    Hmat <- solve(t(wt*bs.old)%*%bs.old + lambda*t(Dmat1)%*%(Dmat1)+
                    kp*t(Dmat)%*%Vmat%*%Dmat)%*%t(wt*bs.old)%*%bs.old
   traceH <- sum(diag(Hmat))
      
    }
    
    ### update beta ####                                                                                                                                                                                                                                                                                                                                                                                                                                          

#       eta.old<- xmats%*%beta.old
       eta.stand <- eta.old
#       q.old <- (eta.stand/atu+1)/2
#       Ut <- pgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )   
#       bs.old <- bs(Ut,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots =c(0,1),intercept=TRUE)
#       
      muhat <- exp(bs.old %*% delta.update)/(1+exp(bs.old %*% delta.update))
      bs.deriv <- splineDesign(knots= c(rep(0,4),knots[c(-1,-length(knots))],rep(1,4)), Ut, ord = 4, derivs=rep(1,length(y0)),outer.ok=TRUE)
      fun.deriv <- bs.deriv%*% delta.update
      den.value<- dgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )/(2*atu)
      du.eta <- size*muhat*(1-muhat)* fun.deriv*den.value
      z.beta <- eta.stand +dd*(y0 - size*muhat)/du.eta
      wt.beta <- as.numeric(size*muhat*(1-muhat)*fun.deriv^2*den.value^2)
       beta.update <- solve(t(wt.beta*xmats)%*%xmats)%*%t(wt.beta*xmats)%*%z.beta
#       cz.beta <- chol(t(xmats)%*%wt.beta%*%xmats)
#       beta.update <- chol2inv(cz.beta)%*%t(xmats)%*%wt.beta%*%z.beta
      if(beta.update[1] >0){
        beta.update <- beta.update/sqrt(sum(beta.update^2))
      }
     if(beta.update[1]<0)
     {
       beta.update <- -beta.update/sqrt(sum(beta.update^2))
     }
    
     diff.total<- sqrt(sum((beta.update - beta.old)^2) + sum((delta.update - delta.old)^2))
     if(diff.total<= tol){break} 
   beta.old <- beta.update
  delta.old <- delta.update
     if(j == MaxIter){print('MaxIter reached without convergence')}
    
  }
  
  
  eta <- xmats%*%beta.update
  q.value <- (eta/atu+1)/2
  Ut <- pgenbeta(q.value,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )
  bs0 <- bs(Ut,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots = c(0,1),intercept=TRUE)
  fitted.values <- exp(bs0%*%delta.update)/(1+exp(bs0%*%delta.update))
  res <- list(eta= eta,est = beta.update,delta=delta.update,fitted.values = fitted.values,
              qv=qv,deg=deg,kp=kp,nknots=nknots,knots=knots,d.value=d.value,atu=atu,
              center=center,sd.value = sd.value,catv=catv,traceH=traceH)
  
  return(res)
}

#### pspline using GCV ###
pspline.gcv3 <- function(y0,xmat,size=1,qv=1,deg = 3,nknots=5,kp=1e6,catv=NULL,
                        monotone=TRUE,beta0,delta0,tol = 1e-8,lamv=seq(5,100,length.out = 20),MaxIter = 100,dd)
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
      
      eta.old <- xmats%*%beta.old
      q.old <- (eta.old/atu+1)/2
      Ut <- pgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )
      
      bs.old <- bs(Ut,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots = c(0,1),intercept=TRUE)
      
      if(!monotone)
      {
        Dmat1 <- matrix(0,ncol = bs.nc, nrow = bs.nc -2)
        Dmat1[cbind(1:(bs.nc-2),1:(bs.nc-2))] <- 1
        Dmat1[cbind(1:(bs.nc-2),2:(bs.nc-1))] <- -2
        Dmat1[cbind(1:(bs.nc-2),3:bs.nc)] <- 1
        
        bs.eta <- bs.old%*%delta.old
        bs.mu <- exp(bs.eta)/(1+exp(bs.eta))
        wt <- as.numeric(size*bs.mu*(1-bs.mu))
        z <- bs.eta+(y0-size*bs.mu)/as.numeric(size*bs.mu*(1-bs.mu)) 
        delta.update <- ginv(t(wt*bs.old)%*%bs.old+lambda*t(Dmat1)%*%(Dmat1))%*%t(wt*bs.old)%*%z
        Hmat <- solve(t(wt*bs.old)%*%bs.old + lambda*t(Dmat1)%*%(Dmat1))%*%t(wt*bs.old)%*%bs.old
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
        wt <- as.numeric(size*bs.mu*(1-bs.mu))
        z <- bs.eta+(y0-size*bs.mu)/as.numeric(size*bs.mu*(1-bs.mu)) 
        Vmat <- diag(as.numeric(Dmat%*%delta.old) <0 ) 
        #         cz <- chol(t(bs.old)%*%wt%*%bs.old + lambda*t(Dmat1)%*%(Dmat1) +kp*t(Dmat)%*%Vmat%*%Dmat)
        #         delta.update <- chol2inv(cz)%*%t(bs.old)%*%wt%*%z
        delta.update <- ginv(t(wt*bs.old)%*%bs.old + lambda*t(Dmat1)%*%(Dmat1) +kp*t(Dmat)%*%Vmat%*%Dmat)%*%t(wt*bs.old)%*%z
        
        Hmat <- solve(t(wt*bs.old)%*%bs.old + lambda*t(Dmat1)%*%(Dmat1)+
                        kp*t(Dmat)%*%Vmat%*%Dmat)%*%t(wt*bs.old)%*%bs.old
        traceH <- sum(diag(Hmat))
        
      }
      
      ### update beta ####                                                                                                                                                                                                                                                                                                                                                                                                                                          
      
      #       eta.old<- xmats%*%beta.old
      eta.stand <- eta.old
      #       q.old <- (eta.stand/atu+1)/2
      #       Ut <- pgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )   
      #       bs.old <- bs(Ut,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots =c(0,1),intercept=TRUE)
      #       
      muhat <- exp(bs.old %*% delta.update)/(1+exp(bs.old %*% delta.update))
      bs.deriv <- splineDesign(knots= c(rep(0,4),knots[c(-1,-length(knots))],rep(1,4)), Ut, ord = 4, derivs=rep(1,length(y0)),outer.ok=TRUE)
      fun.deriv <- bs.deriv%*% delta.update
      den.value<- dgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )/(2*atu)
      du.eta <- size*muhat*(1-muhat)* fun.deriv*den.value
      z.beta <- eta.stand +dd*(y0 - size*muhat)/du.eta
      wt.beta <- as.numeric(size*muhat*(1-muhat)*fun.deriv^2*den.value^2)
      beta.update <- solve(t(wt.beta*xmats)%*%xmats)%*%t(wt.beta*xmats)%*%z.beta
      #       cz.beta <- chol(t(xmats)%*%wt.beta%*%xmats)
      #       beta.update <- chol2inv(cz.beta)%*%t(xmats)%*%wt.beta%*%z.beta
      if(beta.update[1] >0){
        beta.update <- beta.update/sqrt(sum(beta.update^2))
      }
      if(beta.update[1]<0)
      {
        beta.update <- -beta.update/sqrt(sum(beta.update^2))
      }
      
      diff.total<- sqrt(sum((beta.update - beta.old)^2) + sum((delta.update - delta.old)^2))
      if(diff.total<= tol){break} 
      beta.old <- beta.update
      delta.old <- delta.update
      if(j == MaxIter){print('MaxIter reached without convergence')}
      
    }
    
    eta <- xmats%*%beta.update
    fitted.values <- muhat
    gcv <- mean((y0 - size*muhat)^2)/(1-traceH/length(y0))^2
    return(gcv)
  }
  
  
  gcvv <- unlist(lapply(lamv,lam.fun))
  #out.gcv <- optimize(f = lam.fun,interval = lam.interval)
  #lam.value <- out.gcv$minimum
  lam.value <- lamv[which.min(gcvv)]
  
  return(lam.value)
}




##### predict based on pspline link function ######

predict.pspline3 <- function(est.obj,newdata)
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

print(c('psplinelink3','pspline.gcv3','predict.pspline3'))