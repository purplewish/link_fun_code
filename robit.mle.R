

robit.mle <- function(y0,x0,beta0,nu0)
{
  xmat <- cbind(1,x0)
  beta.len <- length(beta0)
  nll <- function(par)
  {
    beta.est <- par[1:2]
    nu.est <- par[3]
    pnu <- pt(q = (xmat%*%beta.est), df = nu.est)
    nll.value <-  - sum(y0*log(pnu) + (1-y0)*log(1 - pnu))
    return(nll.value)
  }
  
  out <- optim(par = c(beta0,nu0),fn = nll,method = 'BFGS',hessian = TRUE)
  est <- out$par
  covm <- solve(out$hessian)
  result <- list(est = est, covm = covm)
  return(result)
}

