###############  one step update in each step (there are two steps) ######
##### use Wang 2015 to estimate monotone functions #######
#### algorithm #### 
#library(actuar)
library(splines)
library(MASS)
library(quadprog)
## qv is value of quantile value of qv ####
# cat is the name of categorical data ###
# when upadting beta, the algorithm does not have alpha gamma

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
psplinelink3<- function(y0,xmat,qv=1,size=1,deg = 3,nknots=10,
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
      bs.mu <- 1/(1+exp(-bs.eta))
      wt <- as.numeric(size*bs.mu*(1-bs.mu))
      index <- wt!=0 & (is.na(wt)==FALSE)
      z <- bs.eta+(y0-size*bs.mu)/wt
      delta.update <- ginv(t(wt[index]*bs.old[index,])%*%bs.old[index,]+lambda*t(Dmat1)%*%(Dmat1))%*%t(wt[index]*bs.old[index,])%*%z[index]
      
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
      bs.mu <- 1/(1+exp(-bs.eta))
      wt <- as.numeric(size*bs.mu*(1-bs.mu))
      z <- y0-size*bs.mu
      
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
    muhat <- 1/(1+exp(-bs.old %*% delta.update))
    bs.deriv <- splineDesign(knots= c(rep(0,4),knots[c(-1,-length(knots))],rep(1,4)), Ut, ord = 4, derivs=rep(1,length(y0)),outer.ok=TRUE)
    fun.deriv <- bs.deriv%*% delta.update
    den.value<- dgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1)/(2*atu)
    du.eta <- size*muhat*(1-muhat)* fun.deriv*den.value
    z.beta <- eta.stand + (y0 - size*muhat)/du.eta
    wt.beta <- as.numeric(size*muhat*(1-muhat)*fun.deriv^2*den.value^2)
    beta.cand <- solve(t(wt.beta*xmats)%*%xmats)%*%t(wt.beta*xmats)%*%z.beta
    
    if(!monotone)
    {
      beta.cand <- beta.cand/sqrt(sum(beta.cand^2))*sign(beta.cand[1])
    }
    else{
      beta.cand <- beta.cand/sqrt(sum(beta.cand^2))
    }
    
    beta.update <- beta.cand
    
    
    
    diff.total<- sqrt(sum((beta.update - beta.old)^2) + sum((delta.update - delta.old)^2))
    
    if(diff.total<= tol)
    { break} 
    beta.old <- beta.update
    delta.old <- delta.update
    
  }
  
  eta <- xmats%*%beta.update
  q.value <- (eta/atu+1)/2
  Ut <- pgenbeta(q.value,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )
  bs0 <- bs(Ut,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots = c(0,1),intercept=TRUE)
  
  if(lambda==0)
  {
    Hmat <- solve(t(wt*bs0[,index])%*%bs0[,index]+ lambda*t(Dmat1)%*%(Dmat1))%*%t(wt*bs0[,index])%*%bs0[,index]
  }else
  {
    Hmat <- ginv(t(wt*bs0)%*%bs.old+ lambda*t(Dmat1)%*%(Dmat1))%*%t(wt*bs0)%*%bs0
  }
  
  traceH <- sum(diag(Hmat))
  
  fitted.values <- 1/(1+exp(-bs0%*%delta.update))
  indicator = 1*(j==MaxIter)
  
  res <- list(eta= eta,est = beta.update,delta=delta.update,fitted.values = fitted.values,
              qv=qv,deg=deg,nknots=nknots,knots=knots,d.value=d.value,atu=atu,
              center=center,sd.value = sd.value,catv=catv,traceH=traceH,message=indicator)
  
  return(res)
}

#### pspline using GCV ###
pspline.gcv3 <- function(y0,xmat,qv=1,size=1,deg = 3,nknots=5,catv=NULL,
                         monotone=TRUE,beta0,delta0,tol = 1e-8,lamv=exp(seq(-5,10,length.out = 50)),MaxIter = 1000)
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
        bs.mu <- 1/(1+exp(-bs.eta))
        wt <- as.numeric(size*bs.mu*(1-bs.mu))
        index <- wt!=0 & (is.na(wt)==FALSE)
        z <- bs.eta+(y0-size*bs.mu)/wt
        if(sum(index)<=1)
        {return(c(0,1))}
        else{
          delta.update <- ginv(t(wt[index]*bs.old[index,])%*%bs.old[index,]+lambda*t(Dmat1)%*%(Dmat1))%*%t(wt[index]*bs.old[index,])%*%z[index]
        }
        
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
        bs.mu <- 1/(1+exp(-bs.eta))
        wt <- as.numeric(size*bs.mu*(1-bs.mu))
        z <- y0-size*bs.mu
        
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
      
      nll.fun <- function(beta.update)
      {
        eta <- xmats%*%beta.update
        q.value <- (eta/atu+1)/2
        Ut <- pgenbeta(q.value,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )
        bs0 <- bs(Ut,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots = c(0,1),intercept=TRUE)
        nllv <- -sum(y0*(bs0%*%delta.update))+sum(size*log(1+exp(bs0%*%delta.update)))
        return(nllv)
      }
      
      eta.stand <- eta.old
      muhat <- 1/(1+exp(-bs.old %*% delta.update))
      bs.deriv <- splineDesign(knots= c(rep(0,4),knots[c(-1,-length(knots))],rep(1,4)), Ut, ord = 4, derivs=rep(1,length(y0)),outer.ok=TRUE)
      fun.deriv <- bs.deriv%*% delta.update
      den.value<- dgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1)/(2*atu)
      du.eta <- size*muhat*(1-muhat)* fun.deriv*den.value
      z.beta <- eta.stand + (y0 - size*muhat)/du.eta
      wt.beta <- as.numeric(size*muhat*(1-muhat)*fun.deriv^2*den.value^2)

      beta.cand <- solve(t(wt.beta*xmats)%*%xmats)%*%t(wt.beta*xmats)%*%z.beta
      
      if(!monotone)
      {
        beta.cand <- beta.cand/sqrt(sum(beta.cand^2))*sign(beta.cand[1])
      }
      else{
        beta.cand <- beta.cand/sqrt(sum(beta.cand^2))
      }
      
      beta.update <- beta.cand

      diff.total<- sqrt(sum((beta.update - beta.old)^2) + sum((delta.update - delta.old)^2))
      
      if(diff.total<= tol)
      { break} 
      beta.old <- beta.update
      delta.old <- delta.update
      
    }
    
    eta <- xmats%*%beta.update
    q.value <- (eta/atu+1)/2
    Ut <- pgenbeta(q.value,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )
    bs0 <- bs(Ut,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots = c(0,1),intercept=TRUE)
    
    if(lambda==0)
    {
      Hmat <- solve(t(wt*bs0[,index])%*%bs0[,index]+ lambda*t(Dmat1)%*%(Dmat1))%*%t(wt*bs0[,index])%*%bs.old[,index]
    }else
    {
      Hmat <- ginv(t(wt*bs0)%*%bs.old+ lambda*t(Dmat1)%*%(Dmat1))%*%t(wt*bs0)%*%bs0
    }
    
    traceH <- sum(diag(Hmat))
    #traceHW <- sum(wt*diag(Hmat))
    
    fitted.values <- 1/(1+exp(-bs0%*%delta.update))
    
    #gcv <- mean((y0 - fitted.values)^2)/(1-traceH/length(y0))^2
    gcv <- -2*sum(y0*log(bs.mu/(1-bs.mu))+log(1-bs.mu))/(1-traceH/length(y0))^2
    # gcv = -sum(y0*log(bs.mu/(1-bs.mu))+log(1-bs.mu)) + traceH/(length(y0)-traceHW)*sum(y0*(y0-bs.mu))
    #gcv <-  -2*sum(y0*log(bs.mu/(1-bs.mu))+log(1-bs.mu)) + 2*traceH/(length(y0)-traceH)*sum((y0-bs.mu)^2/(bs.mu*(1-bs.mu)))
    #gcv <- mean((y0 - fitted.values)^2/wt)/(1-traceH/length(y0))^2
    indicator <- 1*(j==MaxIter)
    return(c(gcv,indicator))
  }
  
  
  gcvv <- lapply(lamv,lam.fun)
  gcvv.mat <- do.call(rbind,gcvv)
  index.lam <- which.min(gcvv.mat[gcvv.mat[,2]==0,1])
  lam.value <- (lamv[gcvv.mat[,2]==0])[index.lam]
  
  return(lam.value)
}



print(c('psplinelink3','pspline.gcv3'))