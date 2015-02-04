
source('gev.mle.R')
library(evd)
library(mgcv)

B <- 1000
phi <- 0.5
ns <- 500
beta0 <- c(1,1)
muhat.logit <- muhat.gev <- muhat.mg <- rep(0,B)
for ( b in 1:B)
{
  set.seed(b)
  x0 <- sort(runif(ns,min = -10,max = -0.1))
  yita0 <- cbind(1,x0)%*%beta0
  prob0 <- 1-pgev(-yita0,loc = 0,scale = 1,shape = phi)
  y0 <- rbinom(ns,size = 1,prob = prob0)
  
  glm.fit.logit <- glm(y0~x0,family = binomial())
  gev.fit <- gev.mle(x0,y = y0,phi = phi,par0 = c(1,1))
  mg.fit<- gam(y0~s(x0,bs = 'cr'),family = binomial())
  
  muhat.logit[b] <- mean((glm.fit.logit$fitted.values-prob0)^2)
  muhat.mg[b] <- mean((mg.fit$fitted.values-prob0)^2) 
  muhat.gev[b] <- mean((gev.fit$fitted.values-prob0)^2) 
}

apply(cbind(muhat.logit,muhat.mg,muhat.gev),2,summary)

plot(x0,prob0,lwd=2,ylim=c(0,1),type = 'l')
lines(x0,gev.fit$fitted.values,col='blue',lwd=2)
lines(x0,glm.fit.logit$fitted.values,col='red',lwd=1)
lines(x0,mg.fit$fitted.values,col='green',lwd=1)







xnew <- runif(100,min=-1,max=-0.1)
yitanew <- cbind(1,xnew)%*%beta0
probnew <- 1-pgev(-yitanew,loc = 0,scale = 1,shape = phi)
ynew <- rbinom(100,size = 1,prob = probnew)

pred.logit <- predict.glm(glm.fit.logit,newdata = data.frame(x0=x0),type='response')
sum(pred.logit[x0>-2]>=0.5)

pred.mg <- predict.gam(mg.fit,newdata = data.frame(x0=x0),type='response')
sum(pred.mg[x0>-2] >= 0.5)
