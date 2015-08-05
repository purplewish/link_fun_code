#### identification solved, no intercept, beta1 =1 #####
### one covariate use wang 2015 algorithm ####

library(splines)
psplinelink50 <- function(y0,x0,deg = 3,lambda,nknots,monotone=FALSE,delta0,tol = 1e-4,boundary,MaxIter=1000)
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
      bs.mu <- exp(bs.eta)/(1+exp(bs.eta))
      wt <- as.numeric(bs.mu*(1-bs.mu))
      z <- bs.eta+(y0-bs.mu)/as.numeric(bs.mu*(1-bs.mu)) 
      delta.update <- solve(t(wt*bs0)%*%bs0 + lambda*t(Dmat1)%*%(Dmat1))%*%t(wt*bs0)%*%z
      diff.value <- sqrt(sum((delta.update-delta.old)^2))
      delta.old <- delta.update
      if(diff.value <= tol){break}
      if(j == MaxIter){print('MaxIter reached without convergence')}      
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
      bs.mu <- exp(bs.eta)/(1+exp(bs.eta))
      wt <- as.numeric(bs.mu*(1-bs.mu))
      z <- y0-bs.mu
    
      index.wt <- wt!=0 & (is.na(wt)==FALSE)
  
      Pmat <- t(wt[index.wt]*bs0[index.wt,])%*%(bs0[index.wt,])+lambda*t(Dmat1)%*%(Dmat1)
      dvec <- t(bs0[index.wt,])%*%(z[index.wt] + (wt[index.wt]*bs0[index.wt,])%*%delta.old)
      
      res.qp <- solve.QP(Dmat=Pmat,dvec=dvec,Amat=t(Dmat),bvec=rep(0,bs.nc-1))
      
      delta.update <- res.qp$solution
      diff.value <- sqrt(sum((delta.update-delta.old)^2))
      delta.old <- delta.update
      if(diff.value <= tol){break} 
      
    }   
    
  }
  
  value <- -sum(y0*log(bs.mu/(1-bs.mu))+log(1-bs.mu))
  indicator <- 1*(j==MaxIter)
  out <- list(fitted.values = bs.mu, eta = bs.eta,value=value,
              delta.est= delta.update,deg=deg,boundary=boundary,knots=knots,message=indicator)
  return(out)
}



####combine together ####
pspline.gcv50 <- function(y0,x0,deg = 3,nknots,monotone=TRUE,delta0,tol = 1e-4,boundary,lamv=seq(5,20,length.out = 10),MaxIter=1000)
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
        bs.mu <- exp(bs.eta)/(1+exp(bs.eta))
        wt <- as.numeric(bs.mu*(1-bs.mu))
        z <- bs.eta+(y0-bs.mu)/as.numeric(bs.mu*(1-bs.mu)) 
        delta.update <- solve(t(wt*bs0)%*%bs0 + lambda*t(Dmat1)%*%(Dmat1))%*%t(wt*bs0)%*%z
        diff.value <- sqrt(sum((delta.update-delta.old)^2))
        delta.old <- delta.update
        if(diff.value <= tol){break}
      }
      
      Hmat <- solve(t(wt*bs0)%*%bs0 + lambda*t(Dmat1)%*%(Dmat1))%*%t(wt*bs0)%*%bs0
      traceH <- sum(diag(Hmat))
      
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
        bs.mu <- exp(bs.eta)/(1+exp(bs.eta))
        wt <- as.numeric(bs.mu*(1-bs.mu))
        z <- y0-bs.mu
        
        index.wt <- wt!=0 & (is.na(wt)==FALSE)
        
        Pmat <- t(wt[index.wt]*bs0[index.wt,])%*%(bs0[index.wt,])+lambda*t(Dmat1)%*%(Dmat1)
        dvec <- t(bs0[index.wt,])%*%(z[index.wt] + (wt[index.wt]*bs0[index.wt,])%*%delta.old)
        
        res.qp <- solve.QP(Dmat=Pmat,dvec=dvec,Amat=t(Dmat),bvec=rep(0,bs.nc-1))
        
        delta.update <- res.qp$solution
        diff.value <- sqrt(sum((delta.update-delta.old)^2))
        delta.old <- delta.update
        if(diff.value <= tol){break} 
      }   
      
      Hmat <- solve(t(wt*bs0)%*%bs0 + lambda*t(Dmat1)%*%(Dmat1))%*%t(wt*bs0)%*%bs0
      traceH <- sum(diag(Hmat))
    }
    
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



predict.pspline50 <- function(est.obj,newdata)
{
  
  deg <- est.obj$deg
  knots <- est.obj$knots
  delta.est <- est.obj$delta.est
  boundary <- est.obj$boundary
  bs.value <- bs(newdata,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots = boundary,intercept=TRUE)
  
  eta.bs <- bs.value%*%delta.est
  prob.est <- exp(eta.bs)/(1+exp(eta.bs))
  return(prob.est)
  
}

print(c('psplinelink50','pspline.gcv50','predict.pspline50'))



