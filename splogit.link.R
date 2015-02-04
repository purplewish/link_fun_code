splogit.link <- function(yita0,r0)
{
  
  if(r0 >0 & r0 <=1)
  {
    prob0 <- exp(yita0)/((1+exp(yita0/r0))^r0)
  }
  
  if(r0 > 1)
  {
    yita.new <- -r0*yita0
    prob0 <- 1-(exp(yita.new)/(1+exp(yita.new)))^(1/r0)
  }
  return(prob0)
}


s1 <- function(para)
{
  exp(para[1])/((1+exp(para[1]/para[2]))^para[2])
}

s2 <- function(para)
{
  yita.new <- -para[2]*para[1]
  prob0 <- 1-(exp(yita.new)/(1+exp(yita.new)))^(1/para[2]) 
  return(prob0)
}

splogit.gr1 <- function(para)
{
  grv1 <- exp(para[1])/(1+exp(para[1]/para[2]))^(para[2]+1)
  ratio.value <- para[1]/para[2]
  grv2.1 <- exp(para[1])/(1+exp(ratio.value))^para[2]
  grv2.2 <- -log(1+exp(ratio.value)) + ratio.value*exp(ratio.value)/(1+exp(ratio.value))
  grv2 <- grv2.1*grv2.2
  grv <- c(grv1,grv2)
  return(grv)
}

splogit.gr2 <- function(para)
{
  prod.value <- -para[2]*para[1]
  grv1 <- exp(-para[1])/(1+exp(prod.value))^{1/para[2]+1}
  grv2.1 <- - exp(-para[1])/(1+exp(prod.value))^(1/para[2])
  grv2.2 <- log(1+exp(prod.value))/para[2]^2 + para[1]*exp(prod.value)/(para[2]*(1+exp(prod.value)))
  grv2 <- grv2.1*grv2.2
  grv <- c(grv1,grv2)
  return(grv)
}


library(numDeriv)
grad(s1,c(3,1))
grad(s2,c(3,1))
splogit.gr1(c(2,1))
splogit.gr2(c(2,1))


s1 <- function(para)
{
  exp(para[1])/(1+exp(para[1]/para[2]))^para[2]
}

s2 <- function(para)
{
  yita.new <- -para[1]*para[2]
  prob0 <- 1/(1+exp(yita.new))^(1/para[2])
  return(prob0)
}
