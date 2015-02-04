## Author: Vivekananda Roy
## Created: Tues, Feb 5, 2013
## Last-Updated: NA
##     Update #: 0
##INPUT
### X  Design matrix
  ### y  Vector of 0s and 1s (if binary data) or the number of successes (if binomial data)
  ### n  Number of trials (needed for binomial data)
###References
###Roy, V. and Kaiser, M. S. (2013) Posterior propriety for Bayesian binomial
###regression modelswith a parametric family of link functions, Statistical
###Methodology, 13: 25--41
###Roy, V. and Hobert, J. P. (2007) Convergence rates and asymptotic
###standard errors for MCMC algorithms for Bayesian probit regression,
### Journal of the Royal Statistical Society, Series B, 69: 607--623
propcheck=function(X,y,n=rep(1,length(y)))
{
  require(boot)

  # Error check input
  n.trials = length(y)
  stopifnot(nrow(X)==n.trials, length(n)==n.trials, all(y<=n))

  # Check rank of X
  if(qr(X)$rank!=ncol(X))
  {
    warning(paste("The posterior is improper because the rank of X is", qr(X)$rank))
    return(FALSE)
  }

  # Set up data structures for binary and binomial observations
  t = rep(1,n.trials)
  if(max(n)==1)
  {
    W  = (1-2*y)*X
  } else
  {
    Xstar=rbind(X,X[y*(y-n)!=0,])
    n.trials=nrow(Xstar)
    t[n==y]=-1
    W=c(t,rep(-1,(n.trials-nrow(X))))*Xstar
  }

  # Check for proper posterior
  tau = rep(1,n.trials)

  res = simplex(a =-tau,
                A1=matrix(1,nr=n.trials,nc=n.trials)-diag(n.trials),
                b1=tau,
                A2=diag(n.trials),
                b2=numeric(n.trials),
                A3=t(W),
                b3=numeric(ncol(X))
               )$soln

  if (min(res)>0)
  {
    return(TRUE)
  } else
  {
    warning("The posterior is not proper because there is no overlap.")
    return(FALSE)
  }
}
