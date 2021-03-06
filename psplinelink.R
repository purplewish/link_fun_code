#### identification solved, no intercept, beta1 =1 #####
library(splines)
psplinelink <- function(y0,x0,size=1,deg = 3,lambda,kp=1e8,nknots,monotone=FALSE,delta0,tol = 1e-4,boundary,MaxIter)
{
  bs.nc <- nknots+deg-1
  beta.old <- 1
  delta.old <- delta0
  
  ##### link function #### 
   
    eta.old <- x0*beta.old
    eta.range <- range(eta.old)
    knots <- seq(eta.range[1],eta.range[2],length.out = nknots)
    bs0 <- bs(eta.old,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots = boundary,intercept=TRUE)
    
    if(monotone==FALSE)
    {
      Dmat1 <- matrix(0,ncol = bs.nc, nrow = bs.nc -2)
      Dmat1[cbind(1:(bs.nc-2),1:(bs.nc-2))] <- 1
      Dmat1[cbind(1:(bs.nc-2),2:(bs.nc-1))] <- -2
      Dmat1[cbind(1:(bs.nc-2),3:bs.nc)] <- 1
      
      for(j in 1:MaxIter)
      {
        bs.eta <- bs0%*%delta.old
        bs.mu <- 1/(1+exp(-bs.eta))
        wt <- as.numeric(bs.mu*(1-bs.mu))
        index.wt <- wt!=0
        z <- bs.eta+(y0-bs.mu)/as.numeric(bs.mu*(1-bs.mu)) 
        if(sum(index.wt)==1)
        {
          delta.update <- ginv(t(wt*bs0)%*%bs0 + lambda*t(Dmat1)%*%(Dmat1))%*%(wt[index.wt]*bs0[index.wt,])%*%z[index.wt]
        }else
        {      delta.update <- ginv(t(wt*bs0)%*%bs0 + lambda*t(Dmat1)%*%(Dmat1))%*%t(wt[index.wt]*bs0[index.wt,])%*%z[index.wt]}
        diff.value <- sqrt(sum((delta.update-delta.old)^2))
        delta.old <- delta.update
        if(diff.value <= tol){break}
      }
      
      if(sum(index.wt)==1)
      {
        Hmat <- ginv(t(wt*bs0)%*%bs0+ lambda*t(Dmat1)%*%(Dmat1))%*%(wt[index.wt]*bs0[index.wt,])%*%bs0[index.wt,]
        
      }else{
        Hmat <- ginv(t(wt*bs0)%*%bs0+ lambda*t(Dmat1)%*%(Dmat1))%*%t(wt[index.wt]*bs0[index.wt,])%*%bs0[index.wt,]
        
      }
      
      
    }
    
    if(monotone==TRUE)
    {
      Dmat <- matrix(0,ncol=bs.nc,nrow=bs.nc-1) ### monotone
      diag(Dmat) <- -1
      Dmat[cbind(1:(bs.nc-1),2:bs.nc)] <- 1
      
      Dmat1 <- matrix(0,ncol = bs.nc, nrow = bs.nc -2)
      Dmat1[cbind(1:(bs.nc-2),1:(bs.nc-2))] <- 1
      Dmat1[cbind(1:(bs.nc-2),2:(bs.nc-1))] <- -2
      Dmat1[cbind(1:(bs.nc-2),3:bs.nc)] <- 1
      
      for(j in 1:MaxIter)
      {
        bs.eta <- bs0%*%delta.old
        bs.mu <- 1/(1+exp(-bs.eta))
        wt <- diag(as.numeric(size*bs.mu*(1-bs.mu)))
        index.wt <- bs.mu*(1-bs.mu)!=0
        z <- bs.eta+(y0-size*bs.mu)/as.numeric(size*bs.mu*(1-bs.mu)) 
        Vmat <- diag(as.numeric(Dmat%*%delta.old) <0 ) 
        delta.update <- ginv(t(bs0)%*%wt%*%bs0 + lambda*t(Dmat1)%*%(Dmat1)+kp*t(Dmat)%*%Vmat%*%Dmat)%*%t(bs0[index.wt,])%*%wt[index.wt,index.wt]%*%z[index.wt]
        diff.value <- sqrt(mean((delta.update-delta.old)^2))
        delta.old <- delta.update
        if(diff.value <= tol){break} 
      }   
      Hmat <- ginv(t(bs0)%*%wt%*%bs0 + lambda*t(Dmat1)%*%(Dmat1)+kp*t(Dmat)%*%Vmat%*%Dmat)%*%t(bs0[index.wt,])%*%wt[index.wt,index.wt]%*%bs0[index.wt,]
     
    }
    
    traceH <- sum(diag(Hmat))
    indicator <- 1*(j==MaxIter)
  value <- -sum(y0*log(bs.mu/(1-bs.mu))+log(1-bs.mu))
  out <- list(fitted.values = bs.mu, eta = bs.eta,value=value,
              delta.est= delta.update,deg=deg,boundary=boundary,knots=knots,traceH=traceH,message=indicator)
  return(out)
}

# AIC.value <- function(est.obj,y0)
# {
#   fitted.value <- est.obj$fitted.value
#   traceH <- est.obj$trace
#   aic.value <- -2*sum(y0*log(fitted.value/(1-fitted.value))+log(1-fitted.value))+2*traceH
#   return(aic.value)
# }


####combine together ####
pspline.gcv <- function(y0,x0,size=1,deg = 3,lam.interval,kp=1e8,nknots,monotone=TRUE,delta0,tol = 1e-4,boundary,lamv=seq(5,20,length.out = 10),MaxIter=1000)
{
  bs.nc <- nknots+deg-1
  lam.fun<- function(lambda)
  {
    beta.old <- 1
    delta.old <- delta0
    
    ##### link function #### 
    
    eta.old <- x0*beta.old
    eta.range <- range(eta.old)
    knots <- seq(eta.range[1],eta.range[2],length.out = nknots)
    bs0 <- bs(eta.old,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots = boundary,intercept=TRUE)
    

    if(monotone==FALSE)
    {
      Dmat1 <- matrix(0,ncol = bs.nc, nrow = bs.nc -2)
      Dmat1[cbind(1:(bs.nc-2),1:(bs.nc-2))] <- 1
      Dmat1[cbind(1:(bs.nc-2),2:(bs.nc-1))] <- -2
      Dmat1[cbind(1:(bs.nc-2),3:bs.nc)] <- 1
      
      for(j in 1:MaxIter)
      {
        bs.eta <- bs0%*%delta.old
        bs.mu <- 1/(1+exp(-bs.eta))
        wt <- as.numeric(bs.mu*(1-bs.mu))
        index.wt <- wt!=0
        z <- bs.eta+(y0-bs.mu)/as.numeric(bs.mu*(1-bs.mu)) 
        if(sum(index.wt)==1)
        {
          delta.update <- ginv(t(wt*bs0)%*%bs0 + lambda*t(Dmat1)%*%(Dmat1))%*%(wt[index.wt]*bs0[index.wt,])%*%z[index.wt]
        }else
        {      delta.update <- ginv(t(wt*bs0)%*%bs0 + lambda*t(Dmat1)%*%(Dmat1))%*%t(wt[index.wt]*bs0[index.wt,])%*%z[index.wt]}
        diff.value <- sqrt(sum((delta.update-delta.old)^2))
        delta.old <- delta.update
        if(diff.value <= tol){break}
      }
      
      if(sum(index.wt)==1)
      {
        Hmat <- ginv(t(wt*bs0)%*%bs0+ lambda*t(Dmat1)%*%(Dmat1))%*%(wt[index.wt]*bs0[index.wt,])%*%bs0[index.wt,]
        
      }else{
        Hmat <- ginv(t(wt*bs0)%*%bs0+ lambda*t(Dmat1)%*%(Dmat1))%*%t(wt[index.wt]*bs0[index.wt,])%*%bs0[index.wt,]
        
      }
      
      
    }
    if(monotone==TRUE)
    {
      Dmat <- matrix(0,ncol=bs.nc,nrow=bs.nc-1) ### monotone
      diag(Dmat) <- -1
      Dmat[cbind(1:(bs.nc-1),2:bs.nc)] <- 1
      
      Dmat1 <- matrix(0,ncol = bs.nc, nrow = bs.nc -2)
      Dmat1[cbind(1:(bs.nc-2),1:(bs.nc-2))] <- 1
      Dmat1[cbind(1:(bs.nc-2),2:(bs.nc-1))] <- -2
      Dmat1[cbind(1:(bs.nc-2),3:bs.nc)] <- 1
      
      for(j in 1:MaxIter)
      {
        bs.eta <- bs0%*%delta.old
        bs.mu <- 1/(1+exp(-bs.eta))
        wt <- diag(as.numeric(size*bs.mu*(1-bs.mu)))
        index.wt <- bs.mu*(1-bs.mu)!=0
        z <- bs.eta+(y0-size*bs.mu)/as.numeric(size*bs.mu*(1-bs.mu)) 
        Vmat <- diag(as.numeric(Dmat%*%delta.old) <0 ) 
        delta.update <- ginv(t(bs0)%*%wt%*%bs0 + lambda*t(Dmat1)%*%(Dmat1)+kp*t(Dmat)%*%Vmat%*%Dmat)%*%t(bs0[index.wt,])%*%wt[index.wt,index.wt]%*%z[index.wt]
        diff.value <- sqrt(mean((delta.update-delta.old)^2))
        delta.old <- delta.update
        if(diff.value <= tol){break} 
      }   
      Hmat <- ginv(t(bs0)%*%wt%*%bs0 + lambda*t(Dmat1)%*%(Dmat1)+kp*t(Dmat)%*%Vmat%*%Dmat)%*%t(bs0[index.wt,])%*%wt[index.wt,index.wt]%*%bs0[index.wt,]
      
    }
    traceH <- sum(diag(Hmat))
    gcv <- mean((y0 - bs.mu)^2)/(1-traceH/length(y0))^2
    indicator <- 1*(j==MaxIter)
    return(c(gcv,indicator))
  }
  
  gcvv <- lapply(lamv,lam.fun)
  gcvv.mat <- do.call(rbind,gcvv)
  index.lam <- which.min(gcvv.mat[gcvv.mat[,2]==0,1])
  lam.value <- (lamv[gcvv.mat[,2]==0])[index.lam]
  
  return(lam.value)
}


predict.pspline <- function(est.obj,newdata)
{
  
  deg <- est.obj$deg
  knots <- est.obj$knots
  delta.est <- est.obj$delta.est
  boundary <- est.obj$boundary
  bs.value <- bs(newdata,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots = boundary,intercept=TRUE)

  eta.bs <- bs.value%*%delta.est
  prob.est <- 1/(1+exp(-eta.bs))
  return(prob.est)
  
}

print(c('psplinelink','pspline.gcv','predict.pspline'))



