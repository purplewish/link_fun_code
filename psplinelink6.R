###############  one step update in each step (there are two steps) ######
##### use nloptr to estimate monotone functions #######
#### algorithm #### 
library(actuar)
library(splines)
library(MASS)
## qv is value of quantile value of qv ####
# cat is the name of categorical data ###

psplinelink6<- function(y0,xmat,qv=1,deg = 3,nknots=10,
                        monotone=TRUE,beta0,delta0,catv=NULL,
                        tol = 1e-8,lambda=20,MaxIter=1000,maxeval=100000)
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
      
      nll.fun<- function(x)
      {
        nllv <-  -sum(y0*(bs.old%*%x))+sum(log(1+exp(bs.old%*%x)))+0.5*lambda*t(x)%*%t(Dmat1)%*%Dmat1%*%x
        gradv <- -colSums(y0*bs.old) + t(bs.old)%*%(exp(bs.old%*%x)/(1+exp(bs.old%*%x)))+lambda*t(Dmat1)%*%Dmat1%*%x
        gradv <- as.numeric(gradv)
        return( list("objective"=nllv, "gradient"= gradv ) )
        
        return(nllv)
      }
      
      eval_g1 <- function(x) 
      {
        return( list("constraints"= -Dmat%*%x,
                     "jacobian"= -Dmat)) 
      }
      
      res <- nloptr( x0=delta0, eval_f=nll.fun, eval_g_ineq = eval_g1, 
                     opts = list("algorithm"="NLOPT_LD_SLSQP", "check_derivatives"=FALSE,maxeval=maxeval))  
      
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
pspline.gcv6 <- function(y0,xmat,qv=1,deg = 3,nknots=5,catv=NULL,
                         monotone=TRUE,beta0,delta0,tol = 1e-8,lamv=exp(seq(-5,10,length.out = 50)),MaxIter = 100,maxeval=100000)
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
        
        nll.fun<- function(x)
        {
          nllv <-  -sum(y0*(bs.old%*%x))+sum(log(1+exp(bs.old%*%x)))+0.5*lambda*t(x)%*%t(Dmat1)%*%Dmat1%*%x
          gradv <- -colSums(y0*bs.old) + t(bs.old)%*%(exp(bs.old%*%x)/(1+exp(bs.old%*%x)))+lambda*t(Dmat1)%*%Dmat1%*%x
          gradv <- as.numeric(gradv)
          return( list("objective"=nllv, "gradient"= gradv ) )
          
          return(nllv)
        }
        
        eval_g1 <- function(x) 
        {
          return( list("constraints"= -Dmat%*%x,
                       "jacobian"= -Dmat)) 
        }
        
        res <- nloptr( x0=delta0, eval_f=nll.fun, eval_g_ineq = eval_g1, 
                       opts = list("algorithm"="NLOPT_LD_SLSQP", "check_derivatives"=FALSE,maxeval=maxeval))  
        
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

predict.pspline6 <- function(est.obj,newdata)
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

print(c('psplinelink6','pspline.gcv6','predict.pspline6'))