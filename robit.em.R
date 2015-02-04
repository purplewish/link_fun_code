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
  
  return(list(beta = beta.new))
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
  yita0 <- as.numeric(xmat%*%beta.new)
  prob0 <- pt(q = yita0,df = nu.new)
  return(list(beta = beta.new, nu = nu.new, fitted.values = prob0))
}

robit.pxem <- function(y0,x0,beta0,nu0,tol=1e-3)
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
    
    out <- optimize(f = nu.function,interval=c(1,10))
    nu.new <- out$minimum
    
    diff.value <- sqrt(sum((beta.old - beta.new)^2)+ (nu.old - nu.new)^2)
    if(diff.value < tol) {break}
    else{beta.old <- beta.new; nu.old <- nu.new}   
  }
  yita0 <- as.numeric(xmat%*%beta.new)
  prob0 <- pt(q = yita0,df = nu.new)
  return(list(beta = beta.new, nu = nu.new, fitted.values = prob0))
  
}

print(c('robit.em','robit.ecmc','robit.pxem'))