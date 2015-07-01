### gev.mle for fixed xi
### gev.mle.xi xi is a parameter
library(nloptr)
library(evd)
gev.mle <- function(x0,y0,xi0,par0=c(0,0))
{
  xmat <- cbind(1,x0)
  nr <- length(y0)
  nll <- function(beta0)
  {
    eta <- xmat%*%beta0
    prob <- 1- pgev(-eta,loc = 0,scale = 1,shape = xi)
    nllv <- sum(y0*(log(prob)-log(1-prob)))+sum(log(1-prob))
    return(-nllv)
  }
  
  gr.gev <- function(beta.est)
  {
    eta.est <- xmat%*%beta.est
    elem <- apply(1- xi*(eta.est),1,function(x){max(0,x)})
    coef.gr <- y0/(1-exp(-elem^(-1/xi)))-1
    gr.beta <- as.numeric(-t(xmat)%*%(coef.gr*elem^(-1/xi-1)))
    return(gr.beta)
  }
  out <- optim(par = c(0.5,0.3),fn = nll,gr=gr.gev,method = 'BFGS',hessian = TRUE)
  
  est <-  out$par
  covm <- solve(out$hessian)
  eta.est <- xmat%*%est
  gr.value <- gr.gev(est)
  fitted.values <- 1- pgev(-eta.est,loc = 0,scale = 1,shape = xi)
  aic <- -2*sum(y0*log(fitted.values/(1-fitted.values))+log(1-fitted.values)) + 2*(length(par0)+1)
  outls <- list(est=est,covm=covm,fitted.values=fitted.values,aic = aic)
  return(outls)
}

##### xi is a parameter to estimate #####
gev.mle.xi <- function(x0,y0,par0=c(0,0,0.5))
{
  
  xmat <- cbind(1,x0)
  nr <- length(y0)
  ncx <- ncol(xmat)
  nll <- function(para)
  {
    beta0 <- para[1:(ncx)]
    xi <- para[ncx+1]
    eta <- xmat%*%beta0
    prob <- 1- pgev(-eta,loc = 0,scale = 1,shape = xi)
    nllv <- sum(y0*(log(prob)-log(1-prob)))+sum(log(1-prob))
    return(-nllv)
  }
  
  gr.gev <- function(para)
  {
    beta.est <- para[1:(ncx)]
    xi <- para[ncx+1]
    eta.est <- xmat%*%beta.est
    elem <- apply(1- xi*(eta.est),1,function(x){max(0,x)})
    coef.gr <- y0/(1-exp(-elem^(-1/xi)))-1
    gr.beta <- -t(xmat)%*%(coef.gr*elem^(-1/xi-1))
    gr.xi <- -sum(coef.gr*(log(elem)/xi^2 + eta.est/(xi*elem))*elem^(-1/xi))
    gr.value <- c(gr.beta,gr.xi)
    return(gr.value)
  }
  out <- optim(par = par0,fn = nll,gr=gr.gev, method = 'BFGS',hessian = TRUE)
  
  est <-  out$par
  covm <- solve(out$hessian)
  eta.est <- xmat%*%(est[1:ncx])
  gr.value <- gr.gev(est)
  fitted.values <- 1- pgev(-eta.est,loc = 0,scale = 1,shape = est[ncx+1])
  
  aic <- -2*sum(y0*log(fitted.values/(1-fitted.values))+log(1-fitted.values)) + 2*(length(par0))
  
  outls <- list(est=est,eta = eta.est,ovm=covm,gr=gr.value,fitted.values=fitted.values,convergence = out$convergence,value=out$value,aic=aic)
  return(outls)
}


###### xi is known, and use nloptr ####
gev.profile <- function(x0,y0,xi,par0=c(0,0),range=c(-1,1),maxeval = 3000)
{
  xmat <- cbind(1,x0)
  nr <- length(y0)
  nll <- function(x,xmat,xi)
  {
    eta <- xmat%*%x
    prob <- 1- pgev(-eta,loc = 0,scale = 1,shape = xi)
    nllv <- -sum(y0*(log(prob)-log(1-prob)))-sum(log(1-prob)) 
    
    elem <- apply(1- xi*(eta),1,function(x){max(0,x)})
    coef.gr <- y0/(1-exp(-elem^(-1/xi)))-1
    gr.beta <- as.numeric(-t(xmat)%*%(coef.gr*elem^(-1/xi-1)))
    return( list("objective"=nllv, 
                 "gradient"= gr.beta ) )
  }
  
  eval_g1 <- function(x,xmat,xi) 
  {
    eta.est <- xmat%*%x
    return( list("constraints"= eta.est*xi - 1,
                 "jacobian"= xi*xmat) )
  }
  
  
  res <- nloptr(x0=par0, eval_f=nll, eval_g_ineq = eval_g1, 
                opts = list("algorithm"="NLOPT_LD_SLSQP", "check_derivatives"=TRUE,maxeval=maxeval), xmat = xmat,xi=xi )  
    
  est <-  res$solution
  eta.est <- xmat%*%est
  gr.gev <- function(beta.est)
  {
    eta.est <- xmat%*%beta.est
    elem <- apply(1- xi*(eta.est),1,function(x){max(0,x)})
    coef.gr <- y0/(1-exp(-elem^(-1/xi)))-1
    gr.beta <- as.numeric(-t(xmat)%*%(coef.gr*elem^(-1/xi-1)))
    return(gr.beta)
  }
  gr.value <- gr.gev(est)
  fitted.values <- 1- pgev(-eta.est,loc = 0,scale = 1,shape = xi)
  outls <- list(est=est,gr=gr.value,fitted.values=fitted.values,value=res$objective,message = res$message)
}

##### use nloptr #### 
gev.mle.new <- function(y0,x0,size=1,par0,maxeval=1000)
{
  xmat <- cbind(1,x0)
  nr <- length(y0) 
  ncx <- ncol(xmat)
  nll <- function(x,xmat)
  {
    eta <- xmat%*%x[1:ncx]
    xi <- x[ncx+1]
    prob <- 1- pgev(-eta,loc = 0,scale = 1,shape = xi)
    nllv <- -sum(y0*(log(prob)-log(1-prob)))-sum(size*log(1-prob)) 
    
    elem <- apply(1- xi*(eta),1,function(x){max(0,x)})
    coef.gr <- y0/(1-exp(-elem^(-1/xi)))-size
    gr.beta <- as.numeric(-t(xmat)%*%(coef.gr*elem^(-1/xi-1)))
    gr.xi <- -sum(coef.gr*(log(elem)/xi^2 + eta/(xi*elem))*elem^(-1/xi))
    gr.value <- c(gr.beta,gr.xi)
    return( list("objective"=nllv, 
                 "gradient"= gr.value ) )
  }
  
  eval_g1 <- function(x, xmat ) 
  {
    eta.est <- xmat%*%x[1:ncx]
    return( list("constraints"= eta.est*x[ncx+1] - 1,
                 "jacobian"= cbind(x[ncx+1]*xmat,eta.est)) )
  }
  
  res <- nloptr( x0=par0, eval_f=nll, eval_g_ineq = eval_g1, 
                 opts = list("algorithm"="NLOPT_LD_SLSQP", "check_derivatives"=TRUE,maxeval=maxeval),xmat = xmat )  
  gr.gev <- function(para)
  {
    beta.est <- para[1:ncx]
    xi <- para[ncx+1]
    eta.est <- xmat%*%beta.est
    elem <- apply(1- xi*(eta.est),1,function(x){max(0,x)})
    coef.gr <- y0/(1-exp(-elem^(-1/xi)))-size
    gr.beta <- -t(xmat)%*%(coef.gr*elem^(-1/xi-1))
    gr.xi <- -sum(coef.gr*(log(elem)/xi^2 + eta.est/(xi*elem))*elem^(-1/xi))
    gr.value <- c(gr.beta,gr.xi)
    return(gr.value)
  }
  
  est <-  res$solution
  eta.est <- xmat%*%est[1:ncx]
  gr.value <- gr.gev(est)
  fitted.values <- 1- pgev(-eta.est,loc = 0,scale = 1,shape = est[ncx+1])
  aic <- -2*sum(y0*log(fitted.values/(1-fitted.values))+log(1-fitted.values)) + 2*(length(par0))
  outls <- list(est=est,eta=eta.est,gr=gr.value,fitted.values=fitted.values,value=res$objective,message = res$message,aic=aic)
  return(outls)
}


predict.gev.new <- function(est.obj,newdata)
{
  ncx <- ncol(newdata)
  est<- est.obj$est
  phiv <- est[ncx+2]
  eta.est <- cbind(1,newdata)%*%est[1:(ncx+1)]
  prob.est <- 1- pgev(-eta.est,loc = 0,scale = 1,shape = phiv)
  return(prob.est)
}

print(c('gev.mle','gev.profile','gev.mle.xi','gev.mle.new','predict.gev.new'))
