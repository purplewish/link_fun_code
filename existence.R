source('code/propcheck.R')
source('code/gev.mle.R')
library(evd)
phi <- -2
n <- 500
x <- runif(n,-5,0.1)
prob <- 1-pgev(q = x,loc = 0,scale = 1,shape = phi)
y <- rbinom(n = n,size = 1,prob = prob)

logl <- function(phi)
{
  value <-sum(-y*(1+phi*x)^(-1/phi) +(1-y)*log(1-exp(-(1+phi*x)^(-1/phi))))
  return(value)
}

x.positv <- x[x>0]
x.neg <- x[x <0] 

phi.seq <- seq(max(-1/x.positv),min(-1/x.neg),by=0.005)

logl.phi <- rep(0,length(phi.seq))
for(i in 1:length(phi.seq))
{
  logl.phi[i] <- logl(phi.seq[i])
}

plot(phi.seq,logl.phi,type='l')
summary(diff(logl.phi))






############ derivative ##########
deriv <- function(phi)
{
  prob <- pgev(q = x,loc = 0,scale = 1,shape = phi)
  value <- sum((log(1+phi*x)/phi^2 - x/(phi*(1+phi*x)))*(y-prob)*log(prob)/(1-prob))
  return(value)
}

x.positv <- x[x>0]
x.neg <- x[x <0]

phi.seq <- seq(max(-1/x.positv),min(-1/x.neg),by=0.005)

deriv.phi <- rep(0,length(phi.seq))
for(i in 1:length(phi.seq))
{
  deriv.phi[i] <- deriv(phi.seq[i])
}

plot(phi.seq,deriv.phi,type='l')

diff(deriv.phi)



###### propcheck####
set.seed(7)
phi <- -2
ns <- 500
beta0 <- c(1,1)
x0 <- sort(runif(ns,min = -1.5,max = 2))
yita0 <- cbind(1,x0)%*%beta0
prob0 <- 1-pgev(-yita0,loc = 0,scale = 1,shape = phi)
y0 <- rbinom(ns,size = 1,prob = prob0)

gev.fit <- gev.mle(x0,y = y0,phi = -2,par0 = c(0,0))
gev.fit1 <- gev.mle.phi(x0,y = y0,par0 = c(0,0,0.5))


source('code/propcheck.R')
propcheck(X = cbind(1,x0),y = y0)
