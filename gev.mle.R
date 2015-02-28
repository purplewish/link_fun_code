### gev.mle for fixed xi
### gev.mle.xi xi is a parameter

gev.mle <- function(x0,y0,xi0,par0=c(0,0))
{
  xmat <- cbind(1,x0)
  nr <- length(y0)
  nll <- function(beta0)
  {
    yita <- xmat%*%beta0
    prob <- 1- pgev(-yita,loc = 0,scale = 1,shape = xi)
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
  yita.est <- xmat%*%est
  gr.value <- gr.gev(est)
  gr.value
  summary(1-yita.est)
  fitted.values <- 1- pgev(-yita.est,loc = 0,scale = 1,shape = xi)
  
  outls <- list(est=est,covm=covm,fitted.values=fitted.values)
  return(outls)
}



###### xi is known, and use nloptr ####
gev.profile <- function(x0,y0,xi,par0=c(0,0),maxeval = 3000)
{
  xmat <- cbind(1,x0)
  nr <- length(y0)
  nll <- function(para,xmat,xi)
  {
    yita <- xmat%*%para
    prob <- 1- pgev(-yita,loc = 0,scale = 1,shape = xi)
    nllv <- -sum(y0*(log(prob)-log(1-prob)))-sum(log(1-prob))
    elem <- apply(1- xi*(yita),1,function(x){max(0,x)})
    coef.gr <- y0/(1-exp(-elem^(-1/xi)))-1
    gr.beta <- as.numeric(-t(xmat)%*%(coef.gr*elem^(-1/xi-1)))
    return( list("objective"=nllv, 
                 "gradient"= gr.beta ) )
  }
  
  eval_g1 <- function(para,xmat,xi) 
  {
    yita.est <- xmat%*%para
    return( list("constraints"= yita.est*xi - 1,
                 "jacobian"= xi*xmat) )
  }
  
  
  res <- nloptr(x0=par0, eval_f=nll, eval_g_ineq = eval_g1, 
                opts = list("algorithm"="NLOPT_LD_SLSQP", "check_derivatives"=TRUE,maxeval=maxeval), xmat = xmat,xi=xi )  
  est <-  res$solution
  yita.est <- xmat%*%est
  gr.gev <- function(beta.est)
  {
    eta.est <- xmat%*%beta.est
    elem <- apply(1- xi*(eta.est),1,function(x){max(0,x)})
    coef.gr <- y0/(1-exp(-elem^(-1/xi)))-1
    gr.beta <- as.numeric(-t(xmat)%*%(coef.gr*elem^(-1/xi-1)))
    return(gr.beta)
  }
  gr.value <- gr.gev(est)
  fitted.values <- 1- pgev(-yita.est,loc = 0,scale = 1,shape = xi)
  outls <- list(est=est,gr=gr.value,fitted.values=fitted.values,value=res$objective,message = res$message)
}

##### xi is a parameter to estimate #####
gev.mle.xi <- function(x0,y0,par0=c(0,0,0.5))
{
  xmat <- cbind(1,x0)
  nr <- length(y0)
  nll <- function(para)
  {
    beta0 <- para[1:2]
    xi <- para[3]
    yita <- xmat%*%beta0
    prob <- 1- pgev(-yita,loc = 0,scale = 1,shape = xi)
    nllv <- sum(y0*(log(prob)-log(1-prob)))+sum(log(1-prob))
    return(-nllv)
  }
  
  gr.gev <- function(para)
  {
    beta.est <- para[1:2]
    xi <- para[3]
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
  yita.est <- xmat%*%(est[1:2])
  gr.value <- gr.gev(est)
  fitted.values <- 1- pgev(-yita.est,loc = 0,scale = 1,shape = est[3])
  
  outls <- list(est=est,eta = yita.est,ovm=covm,gr=gr.value,fitted.values=fitted.values,convergence = out$convergence,value=out$value)
  return(outls)
}

##### use nloptr #### 
gev.mle.new <- function(y0,x0,par0,maxeval=1000)
{
  xmat <- cbind(1,x0)
  nr <- length(y0) 
  nll <- function(x,xmat)
  {
    yita <- xmat%*%x[1:2]
    xi <- x[3]
    prob <- 1- pgev(-yita,loc = 0,scale = 1,shape = xi)
    nllv <- -sum(y0*(log(prob)-log(1-prob)))-sum(log(1-prob)) 
    
    elem <- apply(1- xi*(yita),1,function(x){max(0,x)})
    coef.gr <- y0/(1-exp(-elem^(-1/x[3])))-1
    gr.beta <- as.numeric(-t(xmat)%*%(coef.gr*elem^(-1/xi-1)))
    gr.xi <- -sum(coef.gr*(log(elem)/xi^2 + yita/(xi*elem))*elem^(-1/xi))
    gr.value <- c(gr.beta,gr.xi)
    return( list("objective"=nllv, 
                 "gradient"= gr.value ) )
  }
  
  eval_g1 <- function(x, xmat ) 
  {
    yita.est <- xmat%*%x[1:2]
    return( list("constraints"= yita.est*x[3] - 1,
                 "jacobian"= cbind(x[3]*xmat,yita.est)) )
  }
  
  res <- nloptr( x0=par0, eval_f=nll, eval_g_ineq = eval_g1, 
                 opts = list("algorithm"="NLOPT_LD_SLSQP", "check_derivatives"=TRUE,maxeval=maxeval),
                 xmat = xmat )  
  gr.gev <- function(para)
  {
    beta.est <- para[1:2]
    xi <- para[3]
    eta.est <- xmat%*%beta.est
    elem <- apply(1- xi*(eta.est),1,function(x){max(0,x)})
    coef.gr <- y0/(1-exp(-elem^(-1/xi)))-1
    gr.beta <- -t(xmat)%*%(coef.gr*elem^(-1/xi-1))
    gr.xi <- -sum(coef.gr*(log(elem)/xi^2 + eta.est/(xi*elem))*elem^(-1/xi))
    gr.value <- c(gr.beta,gr.xi)
    return(gr.value)
  }
  
  est <-  res$solution
  yita.est <- xmat%*%est[1:2]
  gr.value <- gr.gev(est)
  fitted.values <- 1- pgev(-yita.est,loc = 0,scale = 1,shape = est[3])
  outls <- list(est=est,eta=yita.est,gr=gr.value,fitted.values=fitted.values,value=res$objective,message = res$message)
  return(outls)
}

print(c('gev.mle','gev.profile','gev.mle.xi','gev.mle.new'))
