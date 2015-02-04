#### profie likelihood for v####
xmat <- cbind(1,x0) 
profile.nu <- function(nu.value, beta0)
{
  nll.function <- function(beta.par)
  {
    eta <- xmat%*%beta.par
    pnu <- pt(q = eta,df = nu.value)
    nll.value <- -sum(y0*log(pnu) + (1-y0)*log(1-pnu))
    return(nll.value)
  }
  out <- optim(par = beta0,fn = nll.function,method = 'BFGS')
  ll.value <- out$value
  return(ll.value)
}

nuv <- seq(0.1,10,0.1)
nllv <- rep(0,length(nuv))
for(i in 1:length(nuv))
{
  nllv[i] <- profile.nu(nuv[i],beta0 = c(1,1))
}

plot(nuv,nllv)
nuv[which.min(nllv)]
