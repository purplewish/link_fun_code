##### splogit #####
## use profile likelihood 
splogit.mle <- function(y0,x0,par0,intervalr)
{
  xmat <- cbind(1,x0)
  splogit.profile <- function(r)
  {
    if(r >0 & r <=1)
    {
      nll.fun <- function(para)
      {
        yita0 <- as.numeric(xmat%*%para)
        nll.value <- sum(y0*yita0 + (1-y0)*log((1+exp(yita0/r))^r-exp(yita0)) -r*log(1+exp(yita0/r)))   
        return(-nll.value)
      }
      
      gr.fun <- function(para)
      {
        yita0 <- as.numeric(xmat%*%para)
        gr1 <- (((1+exp(yita0/r))^(r-1))*exp(yita0/r) - exp(yita0))/((1+exp(yita0/r))^r-exp(yita0))
        gr2 <- exp(yita0/r)/(1+exp(yita0/r))
        gr.value <- -t(xmat)%*%(y0 + (1-y0)*gr1 - gr2)
        return(gr.value)
      }
      
      out <- optim(par = par0,fn = nll.fun,gr = gr.fun,method = 'BFGS')
      
    }
    
    if(r > 1)
    {
      nll.fun <- function(para)
      {
        yita0 <- as.numeric(xmat%*%para)
        yita.new <- -r*yita0
        nll.value <- sum(y0*log((1+exp(yita.new))^(1/r) - exp(-yita0)) - (1-y0)*yita0 - (1/r)*log(1+exp(yita.new)))  
        return(-nll.value)
      }
      
      gr.fun <- function(para)
      {
        yita0 <- as.numeric(xmat%*%para)
        yita.new <- -r*yita0
        gr1.num<- ((1+exp(yita.new))^(1/r-1))* exp(yita.new) - exp(-yita0)
        gr1.den <- (1+exp(yita.new))^(1/r) - exp(-yita0)
        gr1 <- gr1.num/gr1.den
        gr2 <- exp(yita.new)/(1+exp(yita.new))
        gr.value <- as.numeric(t(xmat)%*%(y0*gr1 +(1-y0) - gr2))
        return(gr.value)
      }
      out <- optim(par = par0,fn = nll.fun,gr = gr.fun,method = 'BFGS')   
    }
    nll.value <- out$value
    return(nll.value)
  }
  outr <- optimize(f = splogit.profile,interval = intervalr,tol=1e-4)
  r <- outr$minimum
  
  if(r >0 & r <=1)
  {
    nll.fun <- function(para)
    {
      yita0 <- as.numeric(xmat%*%para)
      nll.value <- sum(y0*yita0 + (1-y0)*log((1+exp(yita0/r))^r-exp(yita0)) -r*log(1+exp(yita0/r)))   
      return(-nll.value)
    }
    
    gr.fun <- function(par)
    {
      yita0 <- as.numeric(xmat%*%par)
      gr1 <- (((1+exp(yita0/r))^(r-1))*exp(yita0/r) - exp(yita0))/((1+exp(yita0/r))^r-exp(yita0))
      gr2 <- exp(yita0/r)/(1+exp(yita0/r))
      gr <- as.numeric(-t(xmat)%*%(y0 + (1-y0)*gr1 - gr2))
      return(gr)
    }
    
    out <- optim(par = par0,fn = nll.fun,gr = gr.fun,method = 'BFGS')
    
  }
  
  if(r > 1)
  {
    nll.fun <- function(par)
    {
      yita0 <- as.numeric(xmat%*%par)
      yita.new <- -r*yita0
      prob0 <- 1-(exp(yita.new)/(1+exp(yita.new)))^(1/r)
      nll.value <- sum(y0*log((1+exp(yita.new))^(1/r) - exp(-yita0)) - (1-y0)*yita0 - (1/r)*log(1+exp(yita.new)))  
      return(-nll.value)
    }
    
    gr.fun <- function(para)
    {
      yita0 <- as.numeric(xmat%*%para)
      yita.new <- -r*yita0
      gr1.num<- ((1+exp(yita.new))^(1/r-1))* exp(yita.new) - exp(-yita0)
      gr1.den <- (1+exp(yita.new))^(1/r) - exp(-yita0)
      gr1 <- gr1.num/gr1.den
      gr2 <- exp(yita.new)/(1+exp(yita.new))
      gr.value <- as.numeric(t(xmat)%*%(y0*gr1 +(1-y0) - gr2))
      return(gr.value)
    }
    out <- optim(par = par0,fn = nll.fun,gr = gr.fun,method = 'BFGS') 
    
  }
  
  fitted.fun <- function(beta.est,r)
  {
    yita0 <- as.numeric(xmat%*%beta.est)
    if(r >0 & r <=1)
    {
      prob0 <- exp(yita0)/((1+exp(yita0/r))^r)
    }
    
    if(r > 1)
    {
      yita.new <- -r*yita0
      prob0 <- 1-(exp(yita.new)/(1+exp(yita.new)))^(1/r)
    }
    return(prob0)
  }
  beta.est <- out$par
  gr.value <- gr.fun(beta.est)
  fitted.values <- fitted.fun(beta.est,r)
  aic <- -2*sum(y0*log(fitted.values/(1-fitted.values))+log(1-fitted.values)) + 2*(length(par0))
  res.eta <- xmat %*% beta.est
  return(list(beta = beta.est, r = r, eta = res.eta, gr=gr.value,convergence = out$convergence,value = out$value,fitted.values=fitted.values,aic=aic))
}








predict.splogit <- function(est.obj,newdata)
{
  beta.est <- est.obj$beta
  r <- est.obj$r

  yita0 <- as.numeric(cbind(1,newdata)%*%beta.est)
  if(r >0 & r <=1)
  {
    prob0 <- exp(yita0)/((1+exp(yita0/r))^r)
  }
  
  if(r > 1)
  {
    yita.new <- -r*yita0
    prob0 <- 1-(exp(yita.new)/(1+exp(yita.new)))^(1/r)
  }

 return(prob0)

}

print(c('splogit.mle','predict.splogit'))