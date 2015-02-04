########### pdf and cdf for different link functions #############
source('burr.fun.R')
library(stats)
library(ggplot2)
library(evd)
library(locfit)
### density function ####

curve(dlogis(x, location = 0, scale = 1),from = -7,to = 7,ylim=c(0,0.45),ylab='density')
curve(dburr(x,shape = 0.5),col='red',add=TRUE)
curve(dgev(-x,loc = 0,scale = 1,shape = -0.5),col='blue',add=TRUE)
curve(dgev(-x,loc = 0,scale = 1,shape = 0.5),col='purple',add=TRUE)

pdf('document/figures/pdf.pdf',width = 10,height = 8)
ggplot(data.frame(x=c(-7, 7)), aes(x)) + 
  stat_function(fun=function(x) dlogis(x,location = 0,scale = 1),aes(color='logit'))+
  stat_function(fun= function(x) dburr(x,shape = 0.5) ,aes(color='burr'))+
  stat_function(fun = function(x) dgev(-x,loc = 0,scale = 1,shape = -0.5),aes(color='gev1'))+
  stat_function(fun = function(x) dgev(-x,loc = 0,scale = 1,shape = 0.5), aes(color='gev2' ))+
  scale_colour_manual(name='link',values=c('logit'='black','burr'="red",'gev1'="blue",'gev2'="green"))
dev.off()

#### cdf ####
pdf('document/figures/cdf.pdf',width = 10,height = 8)
ggplot(data.frame(x=c(-7, 7)), aes(x)) + 
  stat_function(fun=function(x) plogis(x,location = 0,scale = 1),aes(color='logit'))+
  stat_function(fun= function(x) pburr(x,shape = 0.5) ,aes(color='burr'))+
  stat_function(fun = function(x) 1-pgev(-x,loc = 0,scale = 1,shape = -0.5),aes(color='gev1'))+
  stat_function(fun = function(x) 1-pgev(-x,loc = 0,scale = 1,shape = 0.5), aes(color='gev2' ))+
  scale_colour_manual(name='link',values=c('logit'='black','burr'="red",'gev1'="blue",'gev2'="green"))
dev.off()

##### comparison of burr link #####
### different burr ####
pdf('document/figures/burr_cdf.pdf',width = 10,height = 8)
ggplot(data.frame(x=c(-10, 10)), aes(x)) + 
  stat_function(fun=function(x) plogis(x,location = 0,scale = 1),aes(color='logit'))+
  stat_function(fun= function(x) pburr(x,shape = 0.3) ,aes(color='burr1'))+
  stat_function(fun= function(x) pburr(x,shape = 0.5) ,aes(color='burr2'))+
  stat_function(fun= function(x) pburr(x,shape = 0.7) ,aes(color='burr3'))+
  scale_colour_manual(name='link',values=c('logit'='black','burr1'='green','burr2'="red",'burr3'="blue"))
dev.off()

pdf('document/figures/burr_pdf.pdf',width = 10,height = 8)
ggplot(data.frame(x=c(-10, 10)), aes(x)) + 
  stat_function(fun=function(x) dlogis(x,location = 0,scale = 1),aes(color='logit'))+
  stat_function(fun= function(x) dburr(x,shape = 0.3) ,aes(color='burr1'))+
  stat_function(fun= function(x) dburr(x,shape = 0.5) ,aes(color='burr2'))+
  stat_function(fun= function(x) dburr(x,shape = 0.7) ,aes(color='burr3'))+
  scale_colour_manual(name='link',values=c('logit'='black','burr1'='green','burr2'="red",'burr3'="blue"))
dev.off()
###### model #####
source('burr.mle.R')
source('burr.fun.R')
library(mgcv)

seed <- 7
phi <- 0.3
ns <-  500
set.seed(seed+1)
x0 <- sort(rnorm(ns,0,3))
beta0 <- c(1,1)
sim.out <- sim.burr(beta0 = beta0,x0 = x0,phi = phi,seed=seed)
y0 <- sim.out$y0
prob0 <- sim.out$prob0


glm.fit.logit <- glm(y0~x0,family = binomial())
burr.fit <- burr.mle(x0,y = y0,phi = phi,par0 = c(0,0))
mg.fit<- gam(y0~s(x0,k=5,bs = 'cr'),family = binomial())
loc.fit <- locfit(y0~x0,family='binomial',link='logit')

pdf('document/figures/result2000_0.3.pdf')
plot(x0,prob0,lwd=2,ylim=c(0,1),type = 'l')
lines(x0,burr.fit$fitted.values,col='blue',lwd=2)
lines(x0,glm.fit.logit$fitted.values,col='red',lwd=1)
lines(x0,mg.fit$fitted.values,col='green',lwd=1)
lines(x0,fitted(loc.fit),col='purple')
dev.off()


mean((glm.fit.logit$fitted.values-prob0)^2)
mean((mg.fit$fitted.values-prob0)^2)   


###### more than one time simulation ######
ss <- c(200,500,1000,200)
nrep <- 50
kinf <- c(3,4,4,5)
beta0 <- c(1,1)

mse.fun <- function(phi0,x_dist='uniform',ss,nrep,kinf,beta0)
{
  mse <- array(0,dim=c(nrep,3,length(ss)))
  for(i in 1:length(ss))
  {
    for(j in 1:nrep)
    {
      seed <- i*j
      phi <- phi0
      ns <- ss[i]
      set.seed(seed+1)
      if(x_dist=='uniform')
      {
        x0 <- sort(runif(ns,-6,10))
      }

      if(x_dist=='normal')
      {
        x0<- sort(rnorm(ns,0,4))
      }
      sim.out <- sim.burr(beta0 = beta0,x0 = x0,phi = phi,seed=seed)
      y0 <- sim.out$y0
      prob0 <- sim.out$prob0
      
      glm.fit.logit <- glm(y0~x0,family = binomial())
      burr.fit <- burr.mle(x0,y = y0,phi = phi,par0 = c(0,0))
      mg.fit<- gam(y0~s(x0,k=kinf[i],bs = 'cr'),family = binomial())
      
      mse[j,1,i] <-  mean((burr.fit$fitted.values-prob0)^2)
      mse[j,2,i] <- mean((glm.fit.logit$fitted.values-prob0)^2)
      mse[j,3,i]<- mean((mg.fit$fitted.values-prob0)^2)   
    }
  }
  return(mse)
}

mse0.5 <- mse.fun(0.5,ss = ss,nrep = nrep,kinf = kinf,beta0 = beta0)
save(mse0.5,file='document/mse_unif_0.5.RData')
apply(mse0.5,c(2,3),summary)

mse0.3 <- mse.fun(0.3,ss = ss,nrep = nrep,kinf = kinf,beta0 = beta0)
save(mse0.3,file='document/mse_unif_0.3.RData')
apply(mse0.3,c(2,3),summary)
