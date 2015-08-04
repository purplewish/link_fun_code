#### robit link ####
### ECMC algorithm ####


robit.em <- function(y0,x0,beta0,nu0,tol = 1e-3)
{
  beta.old <- beta0
  xmat <- cbind(1,x0)
  repeat
  {
    mu.old <- xmat%*%beta.old
    
    ## E step 
    ## tau
    coef.nu <- (1+2/nu0)^0.5
    prob.mu <- pt(q = - mu.old,df = nu0)
    prob.coef.mu <- pt(q = -coef.nu*mu.old,df = nu0+2)
    tauhat <- (y0 - (2*y0 - 1)*prob.coef.mu)/(y0 - (2*y0 - 1)*prob.mu)
    
    ## z
    zhat <- mu.old + (2*y0 - 1)* dt(x = mu.old,df = nu0)/(y0 - (2*y0 -1)*prob.coef.mu )    
    ## update 
    wt <- diag(as.numeric(tauhat))
    beta.new <- solve(t(xmat) %*% wt %*% xmat) %*% t(xmat) %*% wt %*% zhat
     
    diff.value <- sqrt(sum((beta.old - beta.new)^2))
    if(diff.value < tol) {break}
    else{beta.old <- beta.new}   
  }
  res.eta <- xmat %*% beta.new
  return(list(beta = beta.new,eta=res.eta))
}

robit.ecmc <- function(y0,x0,beta0,nu0,tol = 1e-3)
{
  beta.old <- beta0
  nu.old <- nu0
  xmat <- cbind(1,x0)
  repeat
  {
    mu.old <- xmat%*%beta.old
    
    ## E step 
    ## tau
    coef.nu <- (1+2/nu.old)^0.5
    prob.mu <- pt(q = - mu.old,df = nu.old)
    prob.coef.mu <- pt(q = -coef.nu*mu.old,df = nu.old+2)
    tauhat <- (y0 - (2*y0 - 1)*prob.coef.mu)/(y0 - (2*y0 - 1)*prob.mu)
    
    ## z
    zhat <- mu.old + (2*y0 - 1)* dt(x = mu.old,df = nu.old)/(y0 - (2*y0 -1)*prob.coef.mu )
    
    ## update 
    wt <- diag(as.numeric(tauhat))
    beta.new <- solve(t(xmat) %*% wt %*% xmat) %*% t(xmat) %*% wt %*% zhat
    
    mu.new <- xmat%*%beta.new
    nu.function <- function(nu.par)
    {
      pnu <- pt(q = - mu.new ,df = nu.par)
      nll <-  - sum(y0*log(1-pnu) + (1-y0)*log(pnu))
      return(nll)
    }
    
    out <- optimize(f = nu.function,interval=c(1,10))
    nu.new <- out$minimum
    
    diff.value <- sqrt(sum((beta.old - beta.new)^2)+ (nu.old - nu.new)^2)
    if(diff.value < tol) {break}
    else{beta.old <- beta.new; nu.old <- nu.new}   
  }
  yita0 <- as.numeric(xmat %*% beta.new)
  prob0 <- pt(q = yita0,df = nu.new)
  
  return(list(beta = beta.new, nu = nu.new,eta = yita0,fitted.values = prob0))
}



##### use pexem ####
robit.pxem <- function(y0,x0,beta0,nu0,tol=1e-3,interval.nu)
{
  beta.old <- beta0
  nu.old <- nu0
  nr <- length(y0)
  xmat <- cbind(1,x0)
  repeat
  {
    mu.old <- xmat%*%beta.old
    
    ## E step 
    ## tau
    coef.nu <- (1+2/nu.old)^0.5
    prob.mu <- pt(q = - mu.old,df = nu.old)
    prob.coef.mu <- pt(q = -coef.nu*mu.old,df = nu.old+2)
    tauhat <- (y0 - (2*y0 - 1)*prob.coef.mu)/(y0 - (2*y0 - 1)*prob.mu)
    ## z
    zhat <- mu.old + (2*y0 - 1)* dt(x = mu.old,df = nu.old)/(y0 - (2*y0 -1)*prob.coef.mu )
    wt <- diag(as.numeric(tauhat))
    Stau <- sum(tauhat)
    Stxx <- t(xmat) %*% wt %*% xmat
    Stxz <- t(xmat) %*% wt %*% zhat
    Stzz <- nr*(nu.old+1) - nu.old*sum(tauhat) + sum(tauhat*(2*mu.old*zhat - mu.old^2))
    
    ## update 

    beta.star <- solve(Stxx) %*% Stxz
    alpha.hat <- Stau/nr
    sigma.square <- (Stzz - t(Stxz) %*% solve(Stxx) %*% Stxz)/nr
    beta.new <- as.numeric(alpha.hat/sqrt(sigma.square))*beta.star
    
    mu.new <- xmat%*%beta.new
    nu.function <- function(nu.par)
    {
      pnu <- pt(q = - mu.new ,df = nu.par)
      nll <-  - sum(y0*log(1-pnu) + (1-y0)*log(pnu))
      return(nll)
    }
    
    out <- optimize(f = nu.function,interval=interval.nu)
    nu.new <- out$minimum
    
    diff.value <- sqrt(sum((beta.old - beta.new)^2)+ (nu.old - nu.new)^2)
    if(diff.value < tol) {break}
    else{beta.old <- beta.new; nu.old <- nu.new}   
  }
  yita0 <- as.numeric(xmat%*%beta.new)
  prob0 <- pt(q = yita0,df = nu.new)
  ### aic #### 
  aic <- -2*sum(y0*log(prob0/(1-prob0))+log(1-prob0)) + 2*(length(beta0)+1)
  return(list(beta = beta.new, nu = nu.new, eta = yita0,fitted.values = prob0,aic=aic))
  
}


predict.pxem <- function(est.obj,newdata)
{
  beta.est <- est.obj$beta
  nu.value <- est.obj$nu
  eta.est <- cbind(1,newdata)%*%beta.est
  prob.est <- pt(q = eta.est,df=nu.value)
  return(prob.est)
}


#### for binomial #####
robit.pxem.binom <- function(y0,x0,size=1,beta0,nu0,tol=1e-3,interval.nu)
{

  beta.old <- beta0
  nu.old <- nu0
  nr <- length(y0)
  
  if(length(size)==1){size <- rep(size,nr)}
  
  xmat <- cbind(1,x0)
  repeat
  {
    
    mu.old <- xmat%*%beta.old
    
    ## E step 
    ## tau
    coef.nu <- (1+2/nu.old)^0.5
    prob.mu <- pt(q = - mu.old,df = nu.old)
    prob.coef.mu <- pt(q = -coef.nu*mu.old,df = nu.old+2)
    tauhat <- prob.coef.mu*(size-y0)/prob.mu + y0*(1-prob.coef.mu)/(1-prob.mu)
    tauzhat <- (size-y0)*(mu.old*prob.coef.mu -dt(x = mu.old,df = nu.old))/prob.mu +
     y0*(mu.old*(1-prob.coef.mu)+dt(x = mu.old,df = nu.old))/(1-prob.mu)
    
    ## z
    wt <- as.numeric(tauhat)
    Stau <- sum(tauhat)
    Stxx <- t(wt*xmat) %*% xmat
    Stxz <- t(xmat)%*%tauzhat
    Stzz <- sum(size)*(nu.old+1) - nu.old*Stau + 2*sum(mu.old*tauzhat) - sum(tauhat*mu.old^2)
     
    ## update 
    
    beta.star <- solve(Stxx) %*% Stxz
    alpha.hat <- Stau/sum(size)
    sigma.square <- (Stzz - t(Stxz) %*% solve(Stxx) %*% Stxz)/sum(size)
    beta.new <- as.numeric(alpha.hat/sqrt(sigma.square))*beta.star
    
    mu.new <- xmat%*%beta.new
    nu.function <- function(nu.par)
    {
      pnu <- pt(q = - mu.new ,df = nu.par)
      nll <-  - sum(y0*log(1-pnu) + (size-y0)*log(pnu))
      return(nll)
    }
    
    out <- optimize(f = nu.function,interval=interval.nu)
    nu.new <- out$minimum
    
    diff.value <- sqrt(sum((beta.old - beta.new)^2)+ (nu.old - nu.new)^2)
    if(diff.value < tol) {break}
    else{beta.old <- beta.new; nu.old <- nu.new}   
  }
  yita0 <- as.numeric(xmat%*%beta.new)
  prob0 <- pt(q = yita0,df = nu.new)
  ### aic #### 
  aic <- -2*sum(y0*log(prob0/(1-prob0))+size*log(1-prob0)) + 2*(length(beta0)+1)
  return(list(beta = beta.new, nu = nu.new, eta = yita0,fitted.values = prob0,aic=aic))
  
}
print(c('robit.em','robit.ecmc','robit.pxem','predict.pxem','robit.pxem.binom'))