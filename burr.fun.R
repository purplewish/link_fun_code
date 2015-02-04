#### burr density and probability function #####
pburr <- function(q,shape)
{
  nq <- length(q)
  pburr.value <- rep(0,nq)
  for(i in 1:nq)
  {
    pburr.value[i] <- 1-(1+exp(q[i]))^(-shape)
  }
  return(pburr.value)
}

dburr <- function(x,shape)
{
  nx <- length(x)
  dburr.value <- rep(0,nx)
  for(i in 1:nx)
  {
    dburr.value[i] <- shape*exp(x[i])/(1+exp(x[i]))^(shape+1)
  }
  return(dburr.value)
}


sim.burr <- function(beta0,x0,phi,seed)
{
  set.seed(seed)
  ns <- length(x0)
  xmat <- cbind(1,x0)
  yita0 <- xmat%*%beta0
  prob0 <- 1-(1+exp(yita0))^(-phi)
  y0 <- rbinom(ns,size = 1,prob = prob0)
  return(list(prob0=prob0,y0=y0))
}
