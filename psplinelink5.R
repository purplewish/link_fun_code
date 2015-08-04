###############  one step update in each step (there are two steps) ######
##### use Wang 2015 to estimate monotone functions #######
#### algorithm #### 
library(actuar)
library(splines)
library(MASS)
## qv is value of quantile value of qv ####
# cat is the name of categorical data ###

### if there is false in index, use this function. return the index with numbers, where the flase one is replaced by the close one which can make sure it is increasing.
index.fun <- function(index)
{
  nl <- length(index)
  indexnew <- 1:nl
  indexfalse <- indexnew[index==FALSE]
  indextrue <- indexnew[index]
  min.fun <- function(ind)
  {
    (1:length(indextrue))[which.min(abs(indextrue - ind))]
  }
  indexnew[index==FALSE] <- unlist(lapply(indexfalse,min.fun))
  indexnew[index] <- 1:(length(indextrue))
  return(indexnew)
}
psplinelink5<- function(y0,xmat,qv=1,deg = 3,nknots=10,
                        monotone=TRUE,beta0,delta0,catv=NULL,
                        tol = 1e-8,lambda=20,MaxIter=1000)
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
      index <-  bs.mu > .Machine$double.eps
      wt <- diag(as.numeric(bs.mu[index]*(1-bs.mu[index])))
      z <- bs.eta[index]+(y0[index]-bs.mu[index])/as.numeric(bs.mu[index]*(1-bs.mu[index])) 
      delta.update <- ginv(t(bs.old[index,])%*%wt%*%bs.old[index,]+lambda*t(Dmat1)%*%(Dmat1))%*%t(bs.old[index,])%*%wt%*%z
      
    }
    
    if(monotone)
    {
      Dmat <- matrix(0,ncol=bs.nc,nrow=bs.nc-1)
      diag(Dmat) <- -1
      Dmat[cbind(1:(bs.nc-1),2:bs.nc)] <- 1
      
      Dmat1 <- matrix(0,ncol = bs.nc, nrow = bs.nc-2)
      Dmat1[cbind(1:(bs.nc-2),1:(bs.nc-2))] <- 1
      Dmat1[cbind(1:(bs.nc-2),2:(bs.nc-1))] <- -2
      Dmat1[cbind(1:(bs.nc-2),3:bs.nc)] <- 1
      
      
      bs.eta <- bs.old%*%delta.old
      bs.mu <- exp(bs.eta)/(1+exp(bs.eta))
      wt <- as.numeric(bs.mu*(1-bs.mu))
      z <- y0-bs.mu
      
      if(lambda ==0)
      {
        index.wt <- wt > .Machine$double.eps & (is.na(wt)==FALSE)
        index <- colSums(abs(bs.old[index.wt,]))!=0
        
        np <- sum(index)
        Dmat <- matrix(0,ncol=np,nrow=np-1)
        diag(Dmat) <- -1
        Dmat[cbind(1:(np-1),2:np)] <- 1
        
        Dmat1 <- matrix(0,ncol = np, nrow = np-2)
        Dmat1[cbind(1:(np-2),1:(np-2))] <- 1
        Dmat1[cbind(1:(np-2),2:(np-1))] <- -2
        Dmat1[cbind(1:(np-2),3:np)] <- 1
        
        Pmat <- t(wt[index.wt]*bs.old[index.wt,index])%*%(bs.old[index.wt,index])+lambda*t(Dmat1)%*%(Dmat1)
        dvec <- t(bs.old[index.wt,index])%*%(z[index.wt] + (wt[index.wt]*bs.old[index.wt,index])%*%delta.old[index])
        
        res.qp <- solve.QP(Dmat=Pmat,dvec=dvec,Amat=t(Dmat),bvec=rep(0,np-1))
        
        delta.sol <- res.qp$solution
        if(sum(index)==bs.nc)
        {delta.update=delta.sol}else{
          indexnew <- index.fun(index)
          delta.update <- delta.sol[indexnew]
        }
        
      }
      
      if(lambda!=0)
      {
        index.wt <- wt!=0 & (is.na(wt)==FALSE)
        Pmat <- t(wt[index.wt]*bs.old[index.wt,])%*%(bs.old[index.wt,])+lambda*t(Dmat1)%*%(Dmat1)
        dvec <- t(bs.old[index.wt,])%*%(z[index.wt] + (wt[index.wt]*bs.old[index.wt,])%*%delta.old)
        
        res.qp <- solve.QP(Dmat=Pmat,dvec=dvec,Amat=t(Dmat),bvec=rep(0,bs.nc-1))
        delta.update <- res.qp$solution
      }
      
      
    }


    ### update beta ####                                                                                                                                                                
    eta.stand <- eta.old
     
    muhat <- exp(bs.old %*% delta.update)/(1+exp(bs.old %*% delta.update))
    bs.deriv <- splineDesign(knots= c(rep(0,4),knots[c(-1,-length(knots))],rep(1,4)), Ut, ord = 4, derivs=rep(1,length(y0)),outer.ok=TRUE)
    fun.deriv <- bs.deriv%*% delta.update
    den.value<- dgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1)/(2*atu)
    du.eta <- muhat*(1-muhat)* fun.deriv*den.value
    index.deriv <- du.eta!=0
    z.beta <- eta.stand +(y0 - muhat)/du.eta
    wt.beta <- as.numeric(muhat*(1-muhat)*fun.deriv^2*den.value^2)
    beta.update <- ginv(t(wt.beta[index.deriv]*xmats[index.deriv,])%*%xmats[index.deriv,])%*%t(wt.beta[index.deriv]*xmats[index.deriv,])%*%z.beta[index.deriv]
  
    if(beta.update[1] >0){
      beta.update <- beta.update/sqrt(sum(beta.update^2))
    }
    if(beta.update[1]<0)
    {
      beta.update <- -beta.update/sqrt(sum(beta.update^2))
    }
    
    diff.total<- sqrt(sum((beta.update - beta.old)^2) + sum((delta.update - delta.old)^2))
    if(diff.total<= tol)
      {
      if(lambda==0)
      {
        Hmat <- solve(t(wt*bs.old[,index])%*%bs.old[,index]+ lambda*t(Dmat1)%*%(Dmat1))%*%t(wt*bs.old[,index])%*%bs.old[,index]
      }else
      {
        Hmat <- solve(t(wt*bs.old)%*%bs.old+ lambda*t(Dmat1)%*%(Dmat1))%*%t(wt*bs.old)%*%bs.old
      }

      traceH <- sum(diag(Hmat))
      break} 
    beta.old <- beta.update
    delta.old <- delta.update
    
    if(j == MaxIter)
      {
      
      if(lambda==0)
      {
        Hmat <- solve(t(wt*bs.old[,index])%*%bs.old[,index]+ lambda*t(Dmat1)%*%(Dmat1))%*%t(wt*bs.old[,index])%*%bs.old[,index]
      }else
      {
        Hmat <- solve(t(wt*bs.old)%*%bs.old+ lambda*t(Dmat1)%*%(Dmat1))%*%t(wt*bs.old)%*%bs.old
      }
      traceH <- sum(diag(Hmat))
      
      print('MaxIter reached without convergence')}
  }
  
  
  eta <- xmats%*%beta.update
  q.value <- (eta/atu+1)/2
  Ut <- pgenbeta(q.value,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )
  bs0 <- bs(Ut,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots = c(0,1),intercept=TRUE)
  
  fitted.values <- exp(bs0%*%delta.update)/(1+exp(bs0%*%delta.update))
  res <- list(eta= eta,est = beta.update,delta=delta.update,fitted.values = fitted.values,
              qv=qv,deg=deg,nknots=nknots,knots=knots,d.value=d.value,atu=atu,
              center=center,sd.value = sd.value,catv=catv,traceH=traceH)
  
  return(res)
}

#### pspline using GCV ###
pspline.gcv5 <- function(y0,xmat,qv=1,deg = 3,nknots=5,catv=NULL,
                         monotone=TRUE,beta0,delta0,tol = 1e-8,lamv=exp(seq(-5,10,length.out = 50)),MaxIter = 100)
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
        index <-  bs.mu > .Machine$double.eps
        wt <- diag(as.numeric(bs.mu[index]*(1-bs.mu[index])))
        z <- bs.eta[index]+(y0[index]-bs.mu[index])/as.numeric(bs.mu[index]*(1-bs.mu[index])) 
        delta.update <- ginv(t(bs.old[index,])%*%wt%*%bs.old[index,]+lambda*t(Dmat1)%*%(Dmat1))%*%t(bs.old[index,])%*%wt%*%z
        
      }
      
      if(monotone)
      {
        Dmat <- matrix(0,ncol=bs.nc,nrow=bs.nc-1)
        diag(Dmat) <- -1
        Dmat[cbind(1:(bs.nc-1),2:bs.nc)] <- 1
        
        Dmat1 <- matrix(0,ncol = bs.nc, nrow = bs.nc-2)
        Dmat1[cbind(1:(bs.nc-2),1:(bs.nc-2))] <- 1
        Dmat1[cbind(1:(bs.nc-2),2:(bs.nc-1))] <- -2
        Dmat1[cbind(1:(bs.nc-2),3:bs.nc)] <- 1
        
        
        bs.eta <- bs.old%*%delta.old
        bs.mu <- exp(bs.eta)/(1+exp(bs.eta))
        wt <- as.numeric(bs.mu*(1-bs.mu))
        z <- y0-bs.mu
        
        if(lambda ==0)
        {
          index.wt <- wt > .Machine$double.eps & (is.na(wt)==FALSE)
          index <- colSums(abs(bs.old[index.wt,]))!=0
          
          np <- sum(index)
          Dmat <- matrix(0,ncol=np,nrow=np-1)
          diag(Dmat) <- -1
          Dmat[cbind(1:(np-1),2:np)] <- 1
          
          Dmat1 <- matrix(0,ncol = np, nrow = np-2)
          Dmat1[cbind(1:(np-2),1:(np-2))] <- 1
          Dmat1[cbind(1:(np-2),2:(np-1))] <- -2
          Dmat1[cbind(1:(np-2),3:np)] <- 1
          
          Pmat <- t(wt[index.wt]*bs.old[index.wt,index])%*%(bs.old[index.wt,index])+lambda*t(Dmat1)%*%(Dmat1)
          dvec <- t(bs.old[index.wt,index])%*%(z[index.wt] + (wt[index.wt]*bs.old[index.wt,index])%*%delta.old[index])
          
          res.qp <- solve.QP(Dmat=Pmat,dvec=dvec,Amat=t(Dmat),bvec=rep(0,np-1))
          
          delta.sol <- res.qp$solution
          if(sum(index)==bs.nc)
          {delta.update=delta.sol}else{
            indexnew <- index.fun(index)
            delta.update <- delta.sol[indexnew]
          }
          
        }
        
        if(lambda!=0)
        {
          index.wt <- wt!=0 & (is.na(wt)==FALSE)
          Pmat <- t(wt[index.wt]*bs.old[index.wt,])%*%(bs.old[index.wt,])+lambda*t(Dmat1)%*%(Dmat1)
          dvec <- t(bs.old[index.wt,])%*%(z[index.wt] + (wt[index.wt]*bs.old[index.wt,])%*%delta.old)
          
          res.qp <- solve.QP(Dmat=Pmat,dvec=dvec,Amat=t(Dmat),bvec=rep(0,bs.nc-1))
          delta.update <- res.qp$solution
        }
      }
      
      
      ### update beta ####                                                                                                                                                                
      eta.stand <- eta.old
      
      muhat <- exp(bs.old %*% delta.update)/(1+exp(bs.old %*% delta.update))
      bs.deriv <- splineDesign(knots= c(rep(0,4),knots[c(-1,-length(knots))],rep(1,4)), Ut, ord = 4, derivs=rep(1,length(y0)),outer.ok=TRUE)
      fun.deriv <- bs.deriv%*% delta.update
      den.value<- dgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1)/(2*atu)
      du.eta <- muhat*(1-muhat)* fun.deriv*den.value
      index.deriv <- du.eta!=0
      z.beta <- eta.stand +(y0 - muhat)/du.eta
      wt.beta <- as.numeric(muhat*(1-muhat)*fun.deriv^2*den.value^2)
      beta.update <- ginv(t(wt.beta[index.deriv]*xmats[index.deriv,])%*%xmats[index.deriv,])%*%t(wt.beta[index.deriv]*xmats[index.deriv,])%*%z.beta[index.deriv]
      
      if(beta.update[1] >0){
        beta.update <- beta.update/sqrt(sum(beta.update^2))
      }
      if(beta.update[1]<0)
      {
        beta.update <- -beta.update/sqrt(sum(beta.update^2))
      }
      
      diff.total<- sqrt(sum((beta.update - beta.old)^2) + sum((delta.update - delta.old)^2))
      if(diff.total<= tol)
      {
        if(lambda==0)
        {
          Hmat <- solve(t(wt*bs.old[,index])%*%bs.old[,index]+ lambda*t(Dmat1)%*%(Dmat1))%*%t(wt*bs.old[,index])%*%bs.old[,index]
        }else
        {
          Hmat <- solve(t(wt*bs.old)%*%bs.old+ lambda*t(Dmat1)%*%(Dmat1))%*%t(wt*bs.old)%*%bs.old
        }
        
        traceH <- sum(diag(Hmat))
        break} 
      beta.old <- beta.update
      delta.old <- delta.update
      
      if(j == MaxIter)
      {
        
        if(lambda==0)
        {
          Hmat <- solve(t(wt*bs.old[,index])%*%bs.old[,index]+ lambda*t(Dmat1)%*%(Dmat1))%*%t(wt*bs.old[,index])%*%bs.old[,index]
        }else
        {
          Hmat <- solve(t(wt*bs.old)%*%bs.old+ lambda*t(Dmat1)%*%(Dmat1))%*%t(wt*bs.old)%*%bs.old
        }
        traceH <- sum(diag(Hmat))
        
        print('MaxIter reached without convergence')}
    }
    
    eta <- xmats%*%beta.update
    q.value <- (eta/atu+1)/2
    Ut <- pgenbeta(q.value,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )
    bs0 <- bs(Ut,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots = c(0,1),intercept=TRUE)
    
    fitted.values <- exp(bs0%*%delta.update)/(1+exp(bs0%*%delta.update))
    
    gcv <- mean((y0 - fitted.values)^2)/(1-traceH/length(y0))^2
    return(gcv)
  }
  
  
  gcvv <- unlist(lapply(lamv,lam.fun))
  lam.value <- lamv[which.min(gcvv)]
  
  return(lam.value)
}




##### predict based on pspline link function ######

predict.pspline5 <- function(est.obj,newdata)
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

print(c('psplinelink5','pspline.gcv5','predict.pspline5'))