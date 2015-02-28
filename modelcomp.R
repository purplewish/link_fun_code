### 
source('code/robit.em.R')
source('code/splogit.mle.R')
source('code/splinelink.R')
source('code/gev.mle.R')
library(evd)

#### true model is GEV ####
set.seed(4)
xi <- 1
ns <- 200
beta0 <- c(1,1)
x0 <- sort(runif(ns,min = -10,max = -0.3))
yita0 <- cbind(1,x0)%*%beta0
prob0 <- 1-pgev(-yita0,loc = 0,scale = 1,shape = xi)
y0 <- rbinom(ns,size = 1,prob = prob0)
summary(prob0)
table(y0)

logit.fit <- glm(y0~x0,family = binomial(link='logit'))
beta0 <- logit.fit$coefficients
probit.fit <- glm(y0 ~ x0, family=binomial(link='probit'))
gev.fit <- gev.mle.xi(x = x0,y = y0,par0 = c(0,0,-1))
res1 <- gev.mle.new(y0 = y0,x0 = x0,par0 = c(1,1,1),maxeval = 3000)
robit.fit <- robit.pxem(y0 = y0,x0 = x0,beta0 = beta0,nu0 = 2,tol = 1e-3) 
splogit.fit <- splogit.mle(y0 = y0,x0 = x0,par0 = c(0.1,0),intervalr = c(0,3))

deg <- 3
nknots <- 20
tol <- 1e-8
boundary = range(x0)
bs.nc <- nknots+deg-1
delta0 <- rep(1/bs.nc,bs.nc)


pspline.lam<- spline.aic(y0,x0,deg=deg,nknots = nknots,kp=1e6,lam.interval = c(0,20000) ,delta0 = delta0,boundary = range(x0), monotone = TRUE)

pspline.fit <- splinelink(y0,x0,deg=deg,nknots = nknots,kp=1e6,lambda = pspline.lam,delta0 = delta0,boundary = range(x0), monotone = TRUE)


plot(x0,prob0,type='l')
lines(x0,logit.fit$fitted.values)
lines(x0,probit.fit$fitted.values)
lines(x0,gev.fit$fitted.values)
lines(x0,res1$fitted.values)
lines(x0,robit.fit$fitted.values)
lines(x0,splogit.fit$fitted.values)
lines(x0,pspline.fit$fitted.value)


##### xi =-1 #####
mse.logit <- mse.probit <- mse.gev <- mse.gev.new <- mse.robit <- mse.splogit <- mse.pspline <- ks.logit <- ks.probit <- ks.gev <- ks.robit <- ks.splogit <- ks.pspline <- ks.gev.new <- rep(0,100) 
grv1.n1<- grv2.n1 <- matrix(0,100,3)
splogit.rv.n1 <- matrix(0,100,2)
boundary1.n1 <- boundary2.n1 <- rep(0,100)
for(s in 1:100)
{
  set.seed(s)
  xi <- -1
  ns <- 200
  beta0 <- c(1,1)
  x0 <- sort(runif(ns,min = -1.9,max = 1.5))
  yita0 <- cbind(1,x0)%*%beta0
  prob0 <- 1-pgev(-yita0,loc = 0,scale = 1,shape = xi)
  y0 <- rbinom(ns,size = 1,prob = prob0)
  
  logit.fit <- glm(y0~x0,family = binomial(link='logit'))
  beta0 <- logit.fit$coefficients
  probit.fit <- glm(y0 ~ x0, family=binomial(link='probit'))
  gev.fit <- gev.mle.xi(x = x0,y = y0,par0 = c(0,0,-1))
  res1 <- gev.mle.new(y0 = y0,x0 = x0,par0 = c(0,0,-1),maxeval = 3000)
  robit.fit <- robit.pxem(y0 = y0,x0 = x0,beta0 = c(0,0),nu0 = 2,tol = 1e-3) 
  splogit.fit <- splogit.mle(y0 = y0,x0 = x0,par0 = c(0,0),intervalr = c(0.03,10))
  
  deg <- 3
  nknots <- 20
  tol <- 1e-8
  boundary = range(x0)
  bs.nc <- nknots+deg-1
  delta0 <- rep(1/bs.nc,bs.nc)
  
  pspline.lam<- spline.aic(y0,x0,deg=deg,nknots = nknots,kp=1e6,lam.interval = c(0,20000) ,delta0 = delta0,boundary = range(x0), monotone = TRUE)
  
  pspline.fit <- splinelink(y0,x0,deg=deg,nknots = nknots,kp=1e6,lambda = pspline.lam ,delta0 = delta0,boundary = range(x0), monotone = TRUE)
  
  
  boundary1.n1[s] <-  min(1-gev.fit$est[3]*cbind(1,x0)%*%gev.fit$est[1:2])
  boundary2.n1[s] <- min(1-res1$est[3]*cbind(1,x0)%*%res1$est[1:2])
  
  mse.logit[s] <- mean((logit.fit$fitted.values - prob0)^2)
  mse.probit[s] <- mean((probit.fit$fitted.values - prob0)^2)
  mse.gev[s] <- mean((gev.fit$fitted.values - prob0)^2)
  mse.gev.new[s] <- mean((res1$fitted.values - prob0)^2)
  mse.robit[s] <- mean((robit.fit$fitted.values - prob0)^2)
  mse.splogit[s] <- mean((splogit.fit$fitted.values - prob0)^2)
  mse.pspline[s] <- mean((pspline.fit$fitted.values - prob0)^2)
  
  ks.logit[s] <- ks.test(prob0,logit.fit$fitted.values)$p.value
  ks.probit[s] <- ks.test(prob0,probit.fit$fitted.values)$p.value
  ks.gev[s] <- ks.test(prob0,gev.fit$fitted.values)$p.value
  ks.gev.new[s] <- ks.test(prob0,res1$fitted.values)$p.value
  ks.robit[s] <- ks.test(prob0,robit.fit$fitted.values)$p.value
  ks.splogit[s] <- ks.test(prob0,splogit.fit$fitted.values)$p.value
  ks.pspline[s] <- ks.test(prob0,pspline.fit$fitted.values)$p.value
  
  
  grv1.n1[s,] <- gev.fit$gr
  grv2.n1[s,] <- res1$gr
  splogit.rv.n1[s,] <- splogit.fit$gr
  
  print(s)
  
}

count.fun <- function(pv,cutoff)
{
  sum(pv <=cutoff)
}

mse.mat.n1 <- cbind(mse.logit,mse.probit,mse.gev,mse.gev.new,mse.robit,mse.splogit,mse.pspline)
apply(mse.mat.n1,2,mean)
pmat.n1 <- cbind(ks.logit,ks.probit,ks.gev,ks.gev.new,ks.robit,ks.splogit,ks.pspline)
apply(pmat.n1,2,function(x) count.fun(x,0.05))
mean(mse.logit)
mean(mse.probit)
mean(mse.gev)
mean(mse.robit)
mean(mse.splogit)
mean(mse.pspline)


xi <- -1
ns <- 500
beta0 <- c(1,1)
x0 <- sort(runif(ns,min = -2,max = 1.5))
yita0 <- cbind(1,x0)%*%beta0
prob0 <- 1-pgev(-yita0,loc = 0,scale = 1,shape = xi)
y0 <- rbinom(ns,size = 1,prob = prob0)
source('code/link.compare.R')


##### xi =1 ####
mse.logit <- mse.probit <- mse.gev <- mse.gev.new <- mse.robit <- mse.splogit <- mse.pspline <- ks.logit <- ks.probit <- ks.gev <- ks.robit <- ks.splogit <- ks.pspline <- ks.gev.new <- rep(0,100) 
grv1.p1<- grv2.p1 <- matrix(0,100,3)
splogit.rv.p1 <- matrix(0,100,2)
boundary1.p1 <- boundary2.p1 <- rep(0,100)
for(s in 1:100)
{
  # s <- 84
  set.seed(s)
  xi <- 1
  ns <- 500
  beta0 <- c(1,1)
  x0 <- sort(runif(ns,min = -10,max = -0.3))
  yita0 <- cbind(1,x0)%*%beta0
  prob0 <- 1-pgev(-yita0,loc = 0,scale = 1,shape = xi)
  y0 <- rbinom(ns,size = 1,prob = prob0)
  
  logit.fit <- glm(y0~x0,family = binomial(link='logit'))
  beta0 <- logit.fit$coefficients
  probit.fit <- glm(y0 ~ x0, family=binomial(link='probit'))
  gev.fit <- gev.mle.xi(x = x0,y = y0,par0 = c(1,1,1))
 r <- gev.mle.new(y0 = y0,x0 = x0,par0 = c(0.1,0.1,1),maxeval = 3000)
  robit.fit <- robit.pxem(y0 = y0,x0 = x0,beta0 = c(0,0),nu0 = 2,tol = 1e-3) 
  splogit.fit <- splogit.mle(y0 = y0,x0 = x0,par0 = c(0,0),intervalr = c(0.01,10))
  
  deg <- 3
  nknots <- 20
  tol <- 1e-8
  boundary = range(x0)
  bs.nc <- nknots+deg-1
  delta0 <- rep(1/bs.nc,bs.nc)
  
  pspline.lam<- spline.aic(y0,x0,deg=deg,nknots = nknots,kp=1e6,lam.interval = c(0,20000) ,delta0 = delta0,boundary = range(x0), monotone = TRUE)
  
  pspline.fit <- splinelink(y0,x0,deg=deg,nknots = nknots,kp=1e6,lambda = pspline.lam ,delta0 = delta0,boundary = range(x0), monotone = TRUE)
  
  
  boundary1.p1[s] <-  min(1-gev.fit$est[3]*cbind(1,x0)%*%gev.fit$est[1:2])
  boundary2.p1[s] <- min(1-res1$est[3]*cbind(1,x0)%*%res1$est[1:2])
  
  mse.logit[s] <- mean((logit.fit$fitted.values - prob0)^2)
  mse.probit[s] <- mean((probit.fit$fitted.values - prob0)^2)
  mse.gev[s] <- mean((gev.fit$fitted.values - prob0)^2)
  mse.gev.new[s] <- mean((res1$fitted.values - prob0)^2)
  mse.robit[s] <- mean((robit.fit$fitted.values - prob0)^2)
  mse.splogit[s] <- mean((splogit.fit$fitted.values - prob0)^2)
  mse.pspline[s] <- mean((pspline.fit$fitted.values - prob0)^2)
  
  ks.logit[s] <- ks.test(prob0,logit.fit$fitted.values)$p.value
  ks.probit[s] <- ks.test(prob0,probit.fit$fitted.values)$p.value
  ks.gev[s] <- ks.test(prob0,gev.fit$fitted.values)$p.value
  ks.gev.new[s] <- ks.test(prob0,res1$fitted.values)$p.value
  ks.robit[s] <- ks.test(prob0,robit.fit$fitted.values)$p.value
  ks.splogit[s] <- ks.test(prob0,splogit.fit$fitted.values)$p.value
  ks.pspline[s] <- ks.test(prob0,pspline.fit$fitted.values)$p.value
  
  
  grv1.p1[s,] <- gev.fit$gr
  grv2.p1[s,] <- res1$gr
  splogit.rv.p1[s,] <- splogit.fit$gr
  
  print(s)
  
}

mse.mat.p1 <- cbind(mse.logit,mse.probit,mse.gev,mse.gev.new,mse.robit,mse.splogit,mse.pspline)
apply(mse.mat.p1,2,mean)
pmat.p1 <- cbind(ks.logit,ks.probit,ks.gev,ks.gev.new,ks.robit,ks.splogit,ks.pspline)
apply(pmat.p1,2,function(x) count.fun(x,0.05))


plot(x0,prob0,type='l')
lines(x0,logit.fit$fitted.values)
lines(x0,probit.fit$fitted.values)
lines(x0,gev.fit$fitted.values)
lines(x0,res1$fitted.values)
lines(x0,robit.fit$fitted.values)
lines(x0,splogit.fit$fitted.values)
lines(x0,pspline.fit$fitted.value)



##### link compare gev #####
source('link_fun_code/link.compare.R')
ns0 <- 500
nrep0 <- 100

out.gev1 <- link.compare(model = 'gev',ns = ns0,nrep = nrep0,min.value = -9.5,max.value = -0.3,model.args = list(beta0=c(1,1),xi=1),init.args = list(init = c(0.1,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.gev2 <- link.compare(model = 'gev',ns = ns0,nrep = nrep0,min.value = -1.9,max.value = 1.5,model.args = list(beta0=c(1,1),xi=-1),init.args = list(init = c(0,0),xi0=-1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.logit <- link.compare(model = 'logit',ns = ns0,nrep = nrep0,min.value = -4,max.value = 2,model.args = list(beta0=c(1,1)),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.probit<- link.compare(model = 'probit',ns = ns0,nrep = nrep0,min.value = -2.5,max.value = 1,model.args = list(beta0=c(1,1)),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.robit1<- link.compare(model = 'robit',ns = ns0,nrep = nrep0,min.value = -5,max.value = 4,model.args = list(beta0=c(1,1),nu=1),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.robit2<- link.compare(model = 'robit',ns = ns0,nrep = nrep0,min.value = -5,max.value = 4,model.args = list(beta0=c(1,1),nu=2),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.splogit.05<- link.compare(model = 'splogit',ns = ns0,nrep = nrep0,
                              min.value = -4,max.value = 0,model.args = list(beta0=c(1,1),r=0.5),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.splogit.01<- link.compare(model = 'splogit',ns = ns0,nrep = nrep0,
                              min.value = -4,max.value = 0,model.args = list(beta0=c(1,1),r=0.1),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.splogit.2<- link.compare(model = 'splogit',ns = ns0,nrep = nrep0,
                             min.value = -2,max.value = 2,model.args = list(beta0=c(1,1),r=2),init.args = list(init = c(0,0),xi0=-1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.splogit.5<- link.compare(model = 'splogit',ns = ns0,nrep = nrep0,
                             min.value = -1.2,max.value = 2,model.args = list(beta0=c(1,1),r=5),init.args = list(init = c(0,0),xi0=-1,nu0=2,r0=1,intervalr=c(0.03,10)))


mse.out <-  cbind(out.logit$mse.mat,out.probit$mse.mat,out.robit1$mse.mat,out.robit2$mse.mat, out.gev1$mse.mat,out.gev2$mse.mat,out.splogit.01$mse.mat,out.splogit.05$mse.mat, out.splogit.2$mse.mat,out.splogit.5$mse.mat)

max.out <- cbind(out.logit$max.mat,out.probit$max.mat,out.robit1$max.mat,out.robit2$max.mat, out.gev1$max.mat,out.gev2$max.mat,out.splogit.01$max.mat,out.splogit.05$max.mat, out.splogit.2$max.mat,out.splogit.5$max.mat)

p.out <- cbind(out.logit$pmat,out.probit$pmat,out.robit1$pmat,out.robit2$pmat, out.gev1$pmat,out.gev2$pmat,out.splogit.01$pmat,out.splogit.05$pmat, out.splogit.2$pmat,out.splogit.5$pmat)

gr1.out <- cbind(out.logit$gr1,out.probit$gr1,out.robit1$gr1,out.robit2$gr1, out.gev1$gr1,out.gev2$gr1,out.splogit.01$gr1,out.splogit.05$gr1, out.splogit.2$gr1, out.splogit.5$gr1)

gr2.out <- cbind(out.logit$gr2,out.probit$gr2,out.robit1$gr2,out.robit2$gr2, out.gev1$gr2,out.gev2$gr2,out.splogit.01$gr2,out.splogit.05$gr2, out.splogit.2$gr2, out.splogit.5$gr2)

splogit.rv.mat <- cbind(out.logit$splogit.rv.mat,out.probit$splogit.rv.mat,out.robit1$splogit.rv.mat,out.robit2$splogit.rv.mat, out.gev1$splogit.rv.mat,out.gev2$splogit.rv.mat,out.splogit.01$splogit.rv.mat, out.splogit.05$splogit.rv.matout.splogit.2$splogit.rv.mat, out.splogit.5$splogit.rv.mat)

boundary1.mat <- cbind(out.logit$boundary1,out.probit$boundary1,out.robit1$boundary1,out.robit2$boundary1, out.gev1$boundary1,out.gev2$boundary1,out.splogit.01$boundary1, out.splogit.05$boundary1,out.splogit.2$boundary1, out.splogit.5$boundary1)

boundary2.mat <- cbind(out.logit$boundary2,out.probit$boundary2,out.robit1$boundary2,out.robit2$boundary2, out.gev1$boundary2,out.gev2$boundary2,out.splogit.01$boundary2, out.splogit.05$boundary2,out.splogit.2$boundary2, out.splogit.5$boundary2)

save(mse.out,max.out,p.out,gr1.out,gr2.out,splogit.rv.mat,boundary1.mat,boundary2.mat,file = 'output/output200.RData')


sqrt(matrix(colMeans(mse.out),ncol=7,byrow=TRUE))

  
