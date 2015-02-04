#### mle estimation for Burr family ####
# phi should be given in the link function, link function is known
burr.mle <- function(x,y,phi,par0=c(0,0))
{
  xmat <- cbind(1,x)
  nr <- length(y)
  nll <- function(beta0)
  {
    yita <- xmat%*%beta0
    prob <- 1- (1+exp(yita))^(-phi) 
    nllv <- sum(y*(log(prob)-log(1-prob)))+sum(log(1-prob))
    return(-nllv)
  }
  out <- optim(par = par0,fn = nll,method = 'BFGS',hessian = TRUE)

  est <-  out$par
  covm <- solve(out$hessian)
  yita.est <- xmat%*%est
  fitted.values <- 1-(1+exp(yita.est))^(-phi)
  outls <- list(est=est,covm=covm,fitted.values=fitted.values)
  return(outls)
}


