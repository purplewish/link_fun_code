##### simulation for GEV ####
source('link_fun_code/gev.mle.R')
source('link_fun_code/splinelink.R')
library(nloptr)
library(evd)
library(mgcv)
set.seed(7)
phi <- 0.5
ns <- 500
beta0 <- c(1,1)
x0 <- sort(runif(ns,min = -6,max = -0.1))
yita0 <- cbind(1,x0)%*%beta0
prob0 <- 1-pgev(-yita0,loc = 0,scale = 1,shape = phi)
y0 <- rbinom(ns,size = 1,prob = prob0)


deg <- 3
nknots <- 20
tol <- 1e-8
boundary = range(x0)
bs.nc <- nknots+deg-1
delta0 <- rep(1/bs.nc,bs.nc)

lam.vector <- seq(301,600,length.out =300)
aic.value <- rep(0, length(lam.vector))
for(i in 1:length(lam.vector))
{
  pspline.fit <- link.est(y0,x0,deg=deg,nknots = nknots,kp=1e6,lambda = lam.vector[i] ,delta0 = delta0,boundary = boundary, monotone = TRUE)  
  aic.value[i] <- AIC.value(est.obj = pspline.fit,y0)
}

plot(aic.value)
which.min(aic.value )

pspline.fit <- link.est(y0,x0,deg=deg,nknots = nknots,kp=1e6,lambda = lam.vector[which.min(aic.value)] ,delta0 = delta0,boundary = range(x0), monotone = TRUE)

range(diff(pspline.fit$fitted.value))

sum(diff(pspline.fit$fitted.value) >0 )


# phi should be given in the link function, link function is known
# also should assume if phi is not known
glm.fit.logit <- glm(y0~x0,family = binomial())
gev.fit <- gev.mle(x0,y = y0,phi = phi,par0 = c(0,0))
gev.fit1 <- gev.mle.phi(x0,y = y0,par0 = c(0,0,0.5))
mg.fit<- gam(y0~s(x0,bs = 'ps',m = c(2,2),k = 20),family = binomial())
#loc.fit <- locfit(y0~x0,family='binomial',link='logit')
#loc.fit1 <- locfit(y0~lp(x0),deg=2,alpha=0.6,family = 'binomial',link='logit')
pdf('document/figures/unknown link/gev_500.pdf.pdf')
plot(x0,prob0,lwd=2,ylim=c(0,1),type = 'l')
lines(x0,gev.fit1$fitted.values,col='blue',lwd=2)
lines(x0,glm.fit.logit$fitted.values,col='red',lwd=1)
#lines(x0,mg.fit$fitted.values,col='green',lwd=1)
lines(x0,pspline.fit$fitted.value,col='orange',lwd=2)
dev.off()

mean((glm.fit.logit$fitted.values-prob0)^2)
mean((mg.fit$fitted.values-prob0)^2) 
mean((gev.fit$fitted.values-prob0)^2)


#### logit is the true model ####
set.seed(7)
ns <- 500
beta0 <- c(1,1)
x0 <- sort(runif(ns,min = -6,max = 2))
yita0 <- cbind(1,x0)%*%beta0
prob0 <- exp(yita0)/(1+exp(yita0))
y0 <- rbinom(ns,size = 1,prob = prob0)

deg <- 3
nknots <- 20
tol <- 1e-8
boundary = range(x0)
bs.nc <- nknots+deg-1
delta0 <- rep(1/bs.nc,bs.nc)


lam.vector <- seq(601,900,length.out =300)
aic.value <- rep(0, length(lam.vector))
for(i in 1:length(lam.vector))
{
  pspline.fit <- link.est(y0,x0,deg=deg,nknots = nknots,kp=1e6,lambda = lam.vector[i] ,delta0 = delta0,boundary = boundary, monotone = TRUE)  
  aic.value[i] <- AIC.value(est.obj = pspline.fit,y0)
}

plot(aic.value)
which.min(aic.value )

pspline.fit <- link.est(y0,x0,deg=deg,nknots = nknots,kp=1e6,lambda = lam.vector[which.min(aic.value)] ,delta0 = delta0,boundary = boundary, monotone = TRUE)

range(diff(pspline.fit$fitted.value))

sum(diff(pspline.fit$fitted.value) >0 )


glm.fit.logit <- glm(y0~x0,family = binomial())
gev.fit <- gev.mle(x0,y = y0,phi = phi,par0 = c(0,0))
gev.fit1 <- gev.mle.phi(x0,y = y0,par0 = c(0,0,0.5))
mg.fit<- gam(y0~s(x0,bs = 'ps',m = c(2,2),k = 20),family = binomial())

pdf('document/figures/unknown link/logit_500.pdf')
plot(x0,prob0,lwd=2,ylim=c(0,1),type = 'l')
lines(x0,gev.fit1$fitted.values,col='blue',lwd=2)
lines(x0,glm.fit.logit$fitted.values,col='red',lwd=1)
lines(x0,pspline.fit$fitted.value,col='orange',lwd=2)
dev.off()



#### x is nonlinear relationship 
set.seed(7)
ns <- 500
beta0 <- c(-1,-1,0.1)
x0 <- sort(runif(ns,min = -4 ,max = 4))
yita0 <- cbind(1,x0,x0^2)%*%beta0
prob0 <- exp(yita0)/(1+exp(yita0))
y0 <- rbinom(ns,size = 1,prob = prob0)





######## a new model #####
pdf('document/figures/mixed.pdf')
plot(yita0,prob0,type='l',xlab=expression(eta),ylab='prob')
abline(v=-4)
abline(v=0)
dev.off()




#### plot 
pdf('document/figures/gev.compare.pdf',width = 10,height = 8)
ggplot(data.frame(x=c(-7, 7)), aes(x)) + 
  stat_function(fun = function(x) 1-pgev(-x,loc = 0,scale = 1,shape = -0.5),aes(color='gev1'))+
  stat_function(fun = function(x) 1-pgev(-x,loc = 0,scale = 1,shape = 0.5), aes(color='gev2' ))+
  scale_colour_manual(name='link',values=c('gev1'='black','gev2'="red"))
dev.off()


ggplot(data.frame(x=c(-7, 7)), aes(x)) + 
  stat_function(fun = function(x) 1-pgev(-x,loc = 0,scale = 1,shape = gev.fit.new$est[3]),aes(color='gev1'))+
  stat_function(fun = function(x) plogis(x), aes(color='gev2' ))

##### check gradient value and boundary in gev distribution######
source('code/gev.mle.R')
set.seed(7)

##### xi = -2 ####
xi <- -2
ns <- 1000
beta0 <- c(1,1)
x0 <- sort(runif(ns,min = -1.5,max = 2))
yita0 <- cbind(1,x0)%*%beta0
prob0 <- 1-pgev(-yita0,loc = 0,scale = 1,shape = xi)
y0 <- rbinom(ns,size = 1,prob = prob0)
summary(prob0)
table(y0)

grv1 <- grv2 <- matrix(0,100,3)
boundary1 <- boundary2 <- rep(0,100)
beta_mat1 <- beta_mat2 <- matrix(0,100,3)
value1 <- value2 <- rep(0,100)
for(s in 1:100)
{
  set.seed(s)
  x0 <- sort(runif(ns,min = -1.5,max = 2))
  yita0 <- cbind(1,x0)%*%beta0
  prob0 <- 1-pgev(-yita0,loc = 0,scale = 1,shape = xi)
  y0 <- rbinom(ns,size = 1,prob = prob0)
  res.logit <- glm(y0~x0,family = 'binomial')
  coef(res.logit)
  gev.fit <- gev.mle.xi(x = x0,y = y0,par0 = c(coef(res.logit),-2))
  res1 <- gev.mle.new(y0 = y0,x0 = x0,par0 = c(0.9,0.9,-2.2),maxeval = 3000)
  res1$gr
  grv1[s,] <- gev.fit$gr
  grv2[s,] <- res1$gr
  boundary1[s] <-  min(1-gev.fit$est[3]*cbind(1,x0)%*%gev.fit$est[1:2])
  boundary2[s] <- min(1-res1$est[3]*cbind(1,x0)%*%res1$est[1:2])
  beta_mat1[s,] <- gev.fit$est
  beta_mat2[s,] <- res1$est
  value1[s] <- gev.fit$value
  value2[s] <- res1$value
}

#### xi =-1 ####

set.seed(28)
xi <- -1
ns <- 500
beta0 <- c(1,1)
x0 <- sort(runif(ns,min = -1.9,max = 1.5))
yita0 <- cbind(1,x0)%*%beta0
prob0 <- 1-pgev(-yita0,loc = 0,scale = 1,shape = xi)
y0 <- rbinom(ns,size = 1,prob = prob0)
summary(prob0)
table(y0)

grv1.n1<- grv2.n1 <- matrix(0,100,3)
boundary1.n1 <- boundary2.n1 <- rep(0,100)
beta_mat1.n1 <- beta_mat2.n1 <- matrix(0,100,3)
value1.n1 <- value2.n1 <- rep(0,100)
for(s in 1:100)
{
  set.seed(s)
  x0 <- sort(runif(ns,min = -1.9,max = 1.5))
  yita0 <- cbind(1,x0)%*%beta0
  prob0 <- 1-pgev(-yita0,loc = 0,scale = 1,shape = xi)
  y0 <- rbinom(ns,size = 1,prob = prob0)
  res.logit <- glm(y0~x0,family = 'binomial')
  gev.fit <- gev.mle.xi(x = x0,y = y0,par0 = c(1,1,-1))
  res1 <- gev.mle.new(y0 = y0,x0 = x0,par0 = c(0,0,-1),maxeval = 3000)
  grv1.n1[s,] <- gev.fit$gr
  grv2.n1[s,] <- res1$gr
  boundary1.n1[s] <-  min(1-gev.fit$est[3]*cbind(1,x0)%*%gev.fit$est[1:2])
  boundary2.n1[s] <- min(1-res1$est[3]*cbind(1,x0)%*%res1$est[1:2])
  beta_mat1.n1[s,] <- gev.fit$est
  beta_mat2.n1[s,] <- res1$est
  value1.n1[s] <- gev.fit$value
  value2.n1[s] <- res1$value
}


####xi=2####
set.seed(7)
xi <- 2
ns <- 1000
beta0 <- c(1,1)
x0 <- sort(runif(ns,min = -20,max = -0.55))
yita0 <- cbind(1,x0)%*%beta0
prob0 <- 1-pgev(-yita0,loc = 0,scale = 1,shape = xi)
y0 <- rbinom(ns,size = 1,prob = prob0)
summary(prob0)
table(y0)

gev.fit <- gev.mle.xi(x = x0,y = y0,par0 = c(0,0,-1))
gev.fit$convergence
gev.fit$est
summary(1-gev.fit$est[3]*cbind(1,x0)%*%gev.fit$est[1:2])
gev.fit$est
gev.fit$gr
gev.fit$value

res1 <- gev.mle.new(y0 = y0,x0 = x0,par0 = c(0.5,0.5,1),maxeval = 3000)
summary(1-res1$est[3]*cbind(1,x0)%*%res1$est[1:2])
res1$est
res1$gr
res1$value
res1$message




### xi =1 ####
set.seed(7)
xi <- 1
ns <- 1000
beta0 <- c(1,1)
x0 <- sort(runif(ns,min = -10,max = -0.3))
yita0 <- cbind(1,x0)%*%beta0
prob0 <- 1-pgev(-yita0,loc = 0,scale = 1,shape = xi)
y0 <- rbinom(ns,size = 1,prob = prob0)

grv1 <- grv2 <- matrix(0,100,3)
boundary1 <- boundary2 <- rep(0,100)
beta_mat1 <- beta_mat2 <- matrix(0,100,3)
value1 <- value2 <- rep(0,100)
for(s in 1:100)
{
  set.seed(s)
  x0 <- sort(runif(ns,min = -10,max = -0.3))
  yita0 <- cbind(1,x0)%*%beta0
  prob0 <- 1-pgev(-yita0,loc = 0,scale = 1,shape = xi)
  y0 <- rbinom(ns,size = 1,prob = prob0)
  res.logit <- glm(y0~x0,family = 'binomial')
  coef(res.logit)
  gev.fit <- gev.mle.xi(x = x0,y = y0,par0 = c(coef(res.logit),1))
  res1 <- gev.mle.new(y0 = y0,x0 = x0,par0 = c(coef(res.logit),1),maxeval = 3000)
  grv1[s,] <- gev.fit$gr
  grv2[s,] <- res1$gr
  boundary1[s] <-  min(1-gev.fit$est[3]*cbind(1,x0)%*%gev.fit$est[1:2])
  boundary2[s] <- min(1-res1$est[3]*cbind(1,x0)%*%res1$est[1:2])
  beta_mat1[s,] <- gev.fit$est
  beta_mat2[s,] <- res1$est
  value1[s] <- gev.fit$value
  value2[s] <- res1$value
}

which.min(boundary1)
which.min(boundary2)
which.max(value1-value2)

value2[97]

boundary1[20]
boundary2[20]
value1[20]
value2[20]
beta_mat1[20,]
beta_mat2[20,]
#### function ####



##### profile likelihood ####
source('code/gev.mle.R')
set.seed(84)
xi <- 1
ns <- 500
beta0 <- c(1,1)
x0 <- sort(runif(ns,min = -10,max = -0.3))
yita0 <- cbind(1,x0)%*%beta0
prob0 <- 1-pgev(-yita0,loc = 0,scale = 1,shape = xi)
y0 <- rbinom(ns,size = 1,prob = prob0)


xiv <- seq(-1.5,,0.05)
len.xi <- length(xiv)
grvp <- matrix(0,len.xi,2)
boundaryp  <- rep(0,len.xi)
valuep <- rep(0,len.xi)
for(s in 1:len.xi)
{
  gev.p <- gev.profile(y0 = y0,x0 = x0,xi=xiv[s],par0 = c(0,0),maxeval = 3000)
  grvp[s,] <- gev.p$gr
  boundaryp[s] <-  min(1-xiv[s]*cbind(1,x0)%*%gev.p$est)
  valuep[s] <- gev.p$value
}

plot(xiv,valuep)




################# explore the existence ######
set.seed(28)
xi <- -1
ns <- 500
beta0 <- c(1,1)
x0 <- sort(runif(ns,min = -1.9,max = 1.5))
yita0 <- cbind(1,x0)%*%beta0
prob0 <- 1-pgev(-yita0,loc = 0,scale = 1,shape = xi)
y0 <- rbinom(ns,size = 1,prob = prob0)
summary(prob0)
table(y0)

res.gev1 <- gev.mle.new(y0 = y0,x0 = x0,par0 = c(1,1,-1),maxeval = 3000)
res.gev2 <- gev.mle.new(y0 = y0,x0 = x0,par0 = c(0,0,1),maxeval = 3000)
which.min(1-res.gev2$est[3]*cbind(1,x0)%*%res.gev2$est[1:2])
yita.est <- cbind(1,x0)%*%res.gev$est[1:2]

summary(log(1-pgev(q = -yita.est,loc = 0,scale = 1,shape = res.gev$est[3])))
res.gev$est
res.gev$value
