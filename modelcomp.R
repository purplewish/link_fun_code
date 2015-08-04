### 
source('link_fun_code/robit.em.R')
source('link_fun_code/splogit.mle.R')
source('link_fun_code/psplinelink1.R')
source('link_fun_code/gev.mle.R')
library(evd)

#### true model is GEV ####
set.seed(4)
xi <- 1
ns <- 500
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
gev.fit.new <- gev.mle.new(y0 = y0,x0 = x0,par0 = c(1,1,1),maxeval = 3000)
robit.fit <- robit.pxem(y0 = y0,x0 = x0,beta0 = beta0,nu0 = 2,tol = 1e-3) 
splogit.fit <- splogit.mle(y0 = y0,x0 = x0,par0 = c(0.1,0),intervalr = c(0,3))

deg <- 3
nknots <- 20
tol <- 1e-8
boundary = range(x0)
bs.nc <- nknots+deg-1
delta0 <- rep(1/bs.nc,bs.nc)


pspline.lam<- spline.aic(y0,x0,deg=deg,nknots = nknots,kp=1e6,lam.interval = c(0,40000) ,delta0 = delta0,boundary = range(x0), monotone = TRUE)

pspline.fit <- splinelink(y0,x0,deg=deg,nknots = nknots,kp=1e6,lambda = pspline.lam,delta0 = delta0,boundary = range(x0), monotone = TRUE)


plot(x0,prob0,type='l')
lines(x0,logit.fit$fitted.values)
lines(x0,probit.fit$fitted.values)
lines(x0,gev.fit$fitted.values)
lines(x0,gev.fit.new$fitted.values)
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
  gev.fit.new <- gev.mle.new(y0 = y0,x0 = x0,par0 = c(0,0,-1),maxeval = 3000)
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
  boundary2.n1[s] <- min(1-gev.fit.new$est[3]*cbind(1,x0)%*%gev.fit.new$est[1:2])
  
  mse.logit[s] <- mean((logit.fit$fitted.values - prob0)^2)
  mse.probit[s] <- mean((probit.fit$fitted.values - prob0)^2)
  mse.gev[s] <- mean((gev.fit$fitted.values - prob0)^2)
  mse.gev.new[s] <- mean((gev.fit.new$fitted.values - prob0)^2)
  mse.robit[s] <- mean((robit.fit$fitted.values - prob0)^2)
  mse.splogit[s] <- mean((splogit.fit$fitted.values - prob0)^2)
  mse.pspline[s] <- mean((pspline.fit$fitted.values - prob0)^2)
  
  ks.logit[s] <- ks.test(prob0,logit.fit$fitted.values)$p.value
  ks.probit[s] <- ks.test(prob0,probit.fit$fitted.values)$p.value
  ks.gev[s] <- ks.test(prob0,gev.fit$fitted.values)$p.value
  ks.gev.new[s] <- ks.test(prob0,gev.fit.new$fitted.values)$p.value
  ks.robit[s] <- ks.test(prob0,robit.fit$fitted.values)$p.value
  ks.splogit[s] <- ks.test(prob0,splogit.fit$fitted.values)$p.value
  ks.pspline[s] <- ks.test(prob0,pspline.fit$fitted.values)$p.value
  
  
  grv1.n1[s,] <- gev.fit$gr
  grv2.n1[s,] <- gev.fit.new$gr
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
  gev.fit.new <- gev.mle.new(y0 = y0,x0 = x0,par0 = c(0.1,0.1,1),maxeval = 3000)
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
  boundary2.p1[s] <- min(1-gev.fit.new$est[3]*cbind(1,x0)%*%gev.fit.new$est[1:2])
  
  mse.logit[s] <- mean((logit.fit$fitted.values - prob0)^2)
  mse.probit[s] <- mean((probit.fit$fitted.values - prob0)^2)
  mse.gev[s] <- mean((gev.fit$fitted.values - prob0)^2)
  mse.gev.new[s] <- mean((gev.fit.new$fitted.values - prob0)^2)
  mse.robit[s] <- mean((robit.fit$fitted.values - prob0)^2)
  mse.splogit[s] <- mean((splogit.fit$fitted.values - prob0)^2)
  mse.pspline[s] <- mean((pspline.fit$fitted.values - prob0)^2)
  
  ks.logit[s] <- ks.test(prob0,logit.fit$fitted.values)$p.value
  ks.probit[s] <- ks.test(prob0,probit.fit$fitted.values)$p.value
  ks.gev[s] <- ks.test(prob0,gev.fit$fitted.values)$p.value
  ks.gev.new[s] <- ks.test(prob0,gev.fit.new$fitted.values)$p.value
  ks.robit[s] <- ks.test(prob0,robit.fit$fitted.values)$p.value
  ks.splogit[s] <- ks.test(prob0,splogit.fit$fitted.values)$p.value
  ks.pspline[s] <- ks.test(prob0,pspline.fit$fitted.values)$p.value
  
  
  grv1.p1[s,] <- gev.fit$gr
  grv2.p1[s,] <- gev.fit.new$gr
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
lines(x0,gev.fit.new$fitted.values)
lines(x0,robit.fit$fitted.values)
lines(x0,splogit.fit$fitted.values)
lines(x0,pspline.fit$fitted.value)



##### comparison #####
source('link_fun_code/link.compare.R')
ns0 <- 500
nrep0 <- 100

out.gev1 <- link.compare(model = 'gev',s0=100,ns = ns0,nrep = nrep0,min.value = -10,max.value = -0.3,model.args = list(beta0=c(1,1),xi=1),init.args = list(init = c(0.1,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.gev2 <- link.compare(model = 'gev',ns = ns0,nrep = nrep0,min.value = -1.5,max.value = 1.4,model.args = list(beta0=c(1,1),xi=-1),init.args = list(init = c(0.1,0),xi0=-1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.logit <- link.compare(model = 'logit',ns = ns0,nrep = nrep0,min.value = -4,max.value = 2,model.args = list(beta0=c(1,1)),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.probit<- link.compare(model = 'probit',ns = ns0,nrep = nrep0,min.value = -2.5,max.value = 1,model.args = list(beta0=c(1,1)),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.robit1<- link.compare(model = 'robit',ns = ns0,nrep = nrep0,min.value = -5,max.value = 4,model.args = list(beta0=c(1,1),nu=1),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.robit2<- link.compare(model = 'robit',ns = ns0,nrep = nrep0,min.value = -5,max.value = 4,model.args = list(beta0=c(1,1),nu=2),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.splogit.05<- link.compare(model = 'splogit',ns = ns0,nrep = nrep0,
                              min.value = -4,max.value = 0,model.args = list(beta0=c(1,1),r=0.5),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.splogit.01<- link.compare(model = 'splogit',ns = ns0,nrep = nrep0,
                              min.value = -3.5,max.value = 0,model.args = list(beta0=c(1,1),r=0.1),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.splogit.2<- link.compare(model = 'splogit',ns = ns0,nrep = nrep0,
                             min.value = -2,max.value = 2,model.args = list(beta0=c(1,1),r=2),init.args = list(init = c(0,0),xi0=-1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.splogit.5<- link.compare(model = 'splogit',ns = ns0,nrep = nrep0,
                             min.value = -1.2,max.value = 2,model.args = list(beta0=c(1,1),r=5),init.args = list(init = c(0,0),xi0=-1,nu0=2,r0=1,intervalr=c(0.03,10)))



mse.out <-  cbind(out.logit$mse.mat,out.probit$mse.mat,out.robit1$mse.mat,out.robit2$mse.mat, out.gev1$mse.mat,out.gev2$mse.mat,out.splogit.01$mse.mat,out.splogit.05$mse.mat, out.splogit.2$mse.mat,out.splogit.5$mse.mat)


save(mse.out,max.out,p.out,gr1.out,gr2.out,splogit.rv.mat,boundary1.mat,boundary2.mat,file = 'output/output2000.RData')


source('link_fun_code/tab.fig.fun.R')

load('output/output1000.RData')
res <- tab.fig.fun(mse.out)
res$mse
gg <- res$gp + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=11))
pdf('document/figures/comparison/plot1000_rel.pdf',width = 10,height = 6)
gg
dev.off()


###### comparison of two covariates ####
source('link_fun_code/link.compare.b3.R')
ns0 <- 100
nrep0 <- 10

out.logit1 <- link.compare.b3(model = 'logit',ns = ns0,nrep = nrep0,s0=0,
                            muv = -0.5,model.args = list(beta0=c(0,1,1)),
                            init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.1,10)),lamv=seq(2,200,length.out = 30), spline.control = list(deg = 3,nknots = 10,dd=1),weights.arg=c('equal','both','left','right'))

out.probit<- link.compare.b3(model = 'probit',ns = ns0,nrep = nrep0,muv = -0.5,model.args = list(beta0=c(0,1,1)),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=seq(2,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.robit1<- link.compare.b3(model = 'robit',ns = ns0,nrep = nrep0,muv = -0.5,model.args = list(beta0=c(0,1,1),nu=1),init.args = list(init = c(0,0,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=seq(2,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.robit2<- link.compare.b3(model = 'robit',ns = ns0,nrep = nrep0,muv = -0.5,model.args = list(beta0=c(0,1,1),nu=2),init.args = list(init = c(0,0,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound = 3,lamv=seq(2,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.robit3<- link.compare.b3(model = 'robit',ns = ns0,s0=0,nrep = nrep0,muv = -0.5,model.args = list(beta0=c(0,1,1),nu=0.6),init.args = list(init = c(0,0,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=seq(2,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.gev1 <- link.compare.b3(model = 'gev',ns = ns0,nrep =nrep0,muv = -0.5,s0=0,iter = 1000, model.args = list(beta0=c(0,1,1),xi=1,locv=-1.5),init.args = list(init = c(0,0,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound=3,spline.control = list(deg = 3,nknots = 10,dd=1),lamv=seq(2,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.gev2 <- link.compare.b3(model = 'gev',ns = ns0,nrep = nrep0,muv = -0.5,s0=0, model.args = list(beta0=c(0,1,1),xi=0.5,locv=-1),init.args = list(init = c(0,0,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),spline.control = list(deg = 3,nknots = 10,dd=1),lamv=seq(2,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.gev3 <- link.compare.b3(model = 'gev',ns = ns0,nrep = nrep0,muv = -0.5,model.args = list(beta0=c(0,1,1),xi=-0.5,locv=0),init.args = list(init = c(0,0.1,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),spline.control = list(deg = 3,nknots = 10,dd=1),lamv=seq(2,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.gev4 <- link.compare.b3(model = 'gev',ns = ns0,nrep = nrep0,muv = -0.5,s0=0,model.args = list(beta0=c(0,1,1),xi=-1,locv=1.2),init.args = list(init = c(0,0.1,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),spline.control = list(deg = 3,nknots = 10,dd=1),lamv=seq(2,200,length.out = 30),weights.arg=c('equal','both','left','right'))


out.splogit.02<- link.compare.b3(model = 'splogit',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(beta0=c(0,1,1),r=0.2),init.args = list(init = c(0,0,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.3,10)),lamv=seq(2,200,length.out = 30),weights.arg=c('equal','both','left','right'))


out.splogit.5<- link.compare.b3(model = 'splogit',s0=0,ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(beta0=c(0,1,1),r=5),init.args = list(init = c(0,0,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound = 3,lamv=seq(2,200,length.out = 30),weights.arg=c('equal','both','left','right'))




prmse.out <-  cbind(out.logit$prmse.mat,out.probit$prmse.mat,out.robit3$prmse.mat,out.robit1$prmse.mat,out.robit2$prmse.mat, out.gev1$prmse.mat,out.gev2$prmse.mat,out.gev3$prmse.mat,out.gev4$prmse.mat,out.splogit.02$prmse.mat, out.splogit.5$prmse.mat)

wprmse.out <-  cbind(out.logit$wprmse.mat,out.probit$wprmse.mat,out.robit3$wprmse.mat,out.robit1$wprmse.mat,out.robit2$wprmse.mat, out.gev1$wprmse.mat,out.gev2$wprmse.mat,out.gev3$wprmse.mat,out.gev4$wprmse.mat,out.splogit.02$wprmse.mat, out.splogit.5$wprmse.mat)

# prrmse.out <-  cbind(out.logit$prrmse.mat,out.probit$prrmse.mat,out.robit3$prrmse.mat,out.robit1$prrmse.mat,out.robit2$prrmse.mat, out.gev1$prrmse.mat,out.gev2$prrmse.mat,out.gev3$prrmse.mat,out.gev4$prrmse.mat,out.splogit.02$prrmse.mat, out.splogit.5$prrmse.mat)

col.name <-  c('logit','probit','robit','gev','splogit','gam','pspline')
row.name <- c('logit','probit','robit(0.6)','robit(1)','robit(2)','gev(1)','gev(0.5)','gev(-0.5)',"gev(-1)",'splogit(.2)','splogit(5)')

save(out.logit,out.probit,out.robit1,out.robit2,out.robit3,out.gev1,out.gev2,out.gev3,out.gev4,out.splogit.02,out.splogit.5,file='output/output500_binarynew.RData')


source('link_fun_code/tab.fig.fun.R')
res <- tab.fig.fun(rmse.out,col.name = col.name,row.name = row.name,remove = FALSE)
resp <- tab.fig.fun(prmse.out,col.name=col.name,row.name=row.name,remove=FALSE)
resw <- tab.fig.fun(wprmse.out,col.name=col.name,row.name=row.name,remove=FALSE)
gg1 <- res$gp + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=11))
gg2 <- resp$gp + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=11))
gg3 <- resw$gp + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=11))

pdf('figures/plot2005_rel.pdf',width = 12,height = 5)
grid.arrange(gg1,gg2,ncol=2)
dev.off()

pdf('figures/box_gev_100.pdf',height = 6,width = 10)
boxplot(out.gev1$prmse.mat)
dev.off()
boxplot(out.gev1$rmse.mat)



###########################new two covariates ######################
#-------------------no bound muv =0----------------------###

source('link_fun_code/link.compare.b3.R')
ns0 <- 200
nrep0 <- 100

out.logit <- link.compare.b(model = 'logit',ns = ns0,nrep = nrep0,s0=0,
                            muv = 0,bound=3,model.args = list(beta0=c(0,1,1)),
                            init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.probit<- link.compare.b(model = 'probit',ns = ns0,nrep = nrep0,muv = 0,bound=3,model.args = list(beta0=c(0,1,1)),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.robit1<- link.compare.b(model = 'robit',ns = ns0,nrep = nrep0,muv = 0,bound=3,model.args = list(beta0=c(0,1,1),nu=1),init.args = list(init = c(0,0,0.1),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10)))

out.robit2<- link.compare.b(model = 'robit',ns = ns0,nrep = nrep0,muv = 0,bound=3,
                            model.args = list(beta0=c(0,1,1),nu=2),init.args = list(init = c(0,0,0.1),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10)))

out.robit3<- link.compare.b(model = 'robit',ns = ns0,s0=0,nrep = nrep0,muv = 0,bound=3,model.args = list(beta0=c(0,1,1),nu=0.6),init.args = list(init = c(0,0,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10)))

out.gev1 <- link.compare.b(model = 'gev',ns = ns0,nrep =nrep0,muv = 0,s0=0,iter = 1000, model.args = list(beta0=c(0,1,1),xi=1,locv=-2.5),init.args = list(init = c(0,0.1,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10)),bound=3,spline.control = list(deg = 3,nknots = 10,dd=1))

out.gev2 <- link.compare.b(model = 'gev',ns = ns0,nrep = nrep0,muv = 0,s0=0, model.args = list(beta0=c(0,1,1),xi=0.5,locv=-1.5),init.args = list(init = c(0,0.1,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10)),spline.control = list(deg = 3,nknots = 10,dd=1))

out.gev3 <- link.compare.b(model = 'gev',ns = ns0,nrep = nrep0,muv = 0,model.args = list(beta0=c(0,1,1),xi=-0.5,locv=0),init.args = list(init = c(0,0.1,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10)),spline.control = list(deg = 3,nknots = 10,dd=1))

out.gev4 <- link.compare.b(model = 'gev',ns = ns0,nrep = nrep0,muv = 0,s0=0,model.args = list(beta0=c(0,1,1),xi=-1,locv=0),init.args = list(init = c(0,0.1,0.1),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10)),spline.control = list(deg = 3,nknots = 10,dd=1))


out.splogit.02<- link.compare.b(model = 'splogit',ns = ns0,s0=3,nrep = nrep0,muv=0,model.args = list(beta0=c(0,1,1),r=0.2),init.args = list(init = c(0,0.1,-0.1),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10)))


out.splogit.5<- link.compare.b(model = 'splogit',s0=0,ns = ns0,nrep = nrep0,muv=0,model.args = list(beta0=c(0,1,1),r=5),init.args = list(init = c(0,0,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10)),bound = 3)

rmse.out <-  cbind(out.logit$rmse.mat,out.probit$rmse.mat,out.robit3$rmse.mat,out.robit1$rmse.mat,out.robit2$rmse.mat, out.gev1$rmse.mat,out.gev2$rmse.mat,out.gev3$rmse.mat,out.gev4$rmse.mat,out.splogit.02$rmse.mat, out.splogit.5$rmse.mat)

# rrmse.out <-  cbind(out.logit$rrmse.mat,out.probit$rrmse.mat,out.robit3$rrmse.mat,out.robit1$rrmse.mat,out.robit2$rrmse.mat, out.gev1$rrmse.mat,out.gev2$rrmse.mat,out.gev3$rrmse.mat,out.gev4$rrmse.mat,out.splogit.02$rrmse.mat, out.splogit.5$rrmse.mat)


prmse.out <-  cbind(out.logit$prmse.mat,out.probit$prmse.mat,out.robit3$prmse.mat,out.robit1$prmse.mat,out.robit2$prmse.mat, out.gev1$prmse.mat,out.gev2$prmse.mat,out.gev3$prmse.mat,out.gev4$prmse.mat,out.splogit.02$prmse.mat, out.splogit.5$prmse.mat)

# prrmse.out <-  cbind(out.logit$prrmse.mat,out.probit$prrmse.mat,out.robit3$prrmse.mat,out.robit1$prrmse.mat,out.robit2$prrmse.mat, out.gev1$prrmse.mat,out.gev2$prrmse.mat,out.gev3$prrmse.mat,out.gev4$prrmse.mat,out.splogit.02$prrmse.mat, out.splogit.5$prrmse.mat)

col.name <-  c('logit','probit','robit','gev','splogit','gam','pspline')
row.name <- c('logit','probit','robit(0.6)','robit(1)','robit(2)','gev(1)','gev(0.5)','gev(-0.5)',"gev(-1)",'splogit(.2)','splogit(5)')

save(rmse.out,prmse.out,file='output/output500_binary.RData')
source('link_fun_code/tab.fig.fun.R')
res <- tab.fig.fun(rmse.out,col.name = col.name,row.name = row.name,remove = FALSE)
resp <- tab.fig.fun(prmse.out,col.name=col.name,row.name=row.name,remove=FALSE)
gg1 <- res$gp + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=11))
gg2 <- resp$gp + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=11))

pdf('figures/plot2005_rel.pdf',width = 12,height = 5)
grid.arrange(gg1,gg2,ncol=2)
dev.off()

pdf('figures/box_gev_100.pdf',height = 6,width = 10)
boxplot(out.gev1$prmse.mat)
dev.off()
boxplot(out.gev1$rmse.mat)


#####------------------------------ comparison of nonlinear------------------------------------- ####
source('link_fun_code/link.compare.n.R')
ns0 <- 100
nrep0 <- 100

out.logit <- link.compare.n(model = 'logit',ns = ns0,nrep = nrep0,model.args = list(beta0=c(0,1,1)),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)),bound = 2)

out.probit<- link.compare.n(model = 'probit',ns = ns0,nrep = nrep0,model.args = list(beta0=c(0,1,1)),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.robit1<- link.compare.n(model = 'robit',ns = ns0,nrep = nrep0,model.args = list(beta0=c(0,1,1),nu=1),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.robit2<- link.compare.n(model = 'robit',ns = ns0,nrep = nrep0,model.args = list(beta0=c(0,1,1),nu=2),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)),bound = 2)

out.robit3<- link.compare.n(model = 'robit',ns = ns0,nrep = nrep0,model.args = list(beta0=c(0,1,1),nu=0.6),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.gev1 <- link.compare.n(model = 'gev',ns = ns0,nrep = nrep0,min.value = -10,max.value = -0.4,model.args = list(beta0=c(0,1,1),xi=1),init.args = list(init = c(0,0.1,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10)))

out.gev2 <- link.compare.n(model = 'gev',ns = ns0,nrep = nrep0,min.value = -10,max.value = 0,model.args = list(beta0=c(0,1,1),xi=0.5),init.args = list(init = c(0,0.1,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10)))

out.gev3 <- link.compare.n(model = 'gev',ns = ns0,nrep = nrep0,min.value = -1.5,max.value = 1.5,model.args = list(beta0=c(0,1,1),xi=-0.5),init.args = list(init = c(0,0.1,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10)))

out.gev4 <- link.compare.n(model = 'gev',ns = ns0,s0 = 56,nrep = nrep0,min.value = -1,max.value = 1.5,model.args = list(beta0=c(0,1,1),xi=-1),init.args = list(init = c(0,0.1,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10)),spline.control = list(deg = 3,nknots = 10,beta0s=c(1,-1)))


#out.splogit.05<- link.compare.b(model = 'splogit',ns = ns0,nrep = nrep0, model.args = list(beta0=c(0,1,1),r=0.5),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)))

out.splogit.02<- link.compare.n(model = 'splogit',ns = ns0,nrep = nrep0,
                                model.args = list(beta0=c(0,1,1),r=0.2),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)),bound=2.3)

#out.splogit.2<- link.compare.b(model = 'splogit',ns = ns0,nrep = nrep0,
#model.args = list(beta0=c(0,1,1),r=2),init.args = list(init = c(0,0,0),xi0=-1,nu0=2,r0=1,intervalr=c(0.03,10)),bound=2.5)

out.splogit.5<- link.compare.n(model = 'splogit',s0=0,ns = ns0,nrep = nrep0,,model.args = list(beta0=c(0,1,1),r=5),init.args = list(init = c(0,0,0),xi0=-1,nu0=2,r0=1,intervalr=c(0.03,10)),bound=2)

mse.out <-  cbind(out.logit$mse.mat,out.robit3$mse.mat,out.robit1$mse.mat,out.robit2$mse.mat, out.gev1$mse.mat,out.gev2$mse.mat,out.gev3$mse.mat)
col.name <-  c('logit','probit','robit','gev','splogit','gam','pspline')
row.name <- c('logit','robit(0.6)','robit(1)','robit(2)','gev(1)','gev(0.5)','gev(-0.5)')

save(mse.out,file='output/outputall100n.RData')
source('link_fun_code/tab.fig.fun.R')
res <- tab.fig.fun(mse.out,col.name = col.name,row.name = row.name,remove = FALSE)
gg <- res$gp + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=11))
pdf('figures/plot100n_rel.pdf',width = 10,height = 6)
gg
dev.off()




#### ------------nonlinear of x1--------#######
source('link_fun_code/link.compare.n1.R')
ns0 <- 100
nrep0 <- 100

out.logit <- link.compare.n1(model = 'logit',ns = ns0,s0=0,nrep = nrep0,muv =-0.5,init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)),bound = 3,lamv = seq(1,50,length.out = 20))

out.probit<- link.compare.n1(model = 'probit',ns = ns0,nrep = nrep0,muv =-0.5,init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)),bound=3,lamv = seq(1,50,length.out = 20))

out.robit1<- link.compare.n1(model = 'robit',ns = ns0,nrep = nrep0,muv =-0.5,model.args = list(nu=1),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)),bound=3,lamv = seq(1,50,length.out = 20))

out.robit2<- link.compare.n1(model = 'robit',ns = ns0,nrep = nrep0,muv =-0.5,model.args = list(nu=2),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)),bound = 3,lamv = seq(1,50,length.out = 20))

out.robit3<- link.compare.n1(model = 'robit',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(nu=0.6),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)),bound=3,lamv = seq(1,50,length.out = 20))

 out.gev1 <- link.compare.n1(model = 'gev',ns = ns0,s0=0,nrep = nrep0,muv=-0.5,model.args = list(xi=1,locv=-2),init.args = list(init = c(0,0),xi0=0.6,nu0=2,r0=1,intervalr=c(0.03,10)),lamv = seq(1,50,length.out = 20))

out.gev2 <- link.compare.n1(model = 'gev',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(xi=0.5,locv=-1),init.args = list(init = c(0,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10)),lamv = seq(1,50,length.out = 20))

out.gev3 <- link.compare.n1(model = 'gev',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(xi=-0.5,locv=0),init.args = list(init = c(0,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10)),lamv = seq(1,50,length.out = 20))

out.gev4 <- link.compare.n1(model = 'gev',ns = ns0,s0 = 0,muv=-0.5,nrep = nrep0,model.args = list(xi=-1,locv=2),init.args = list(init = c(0,0.2),xi0=-0.5,nu0=1,r0=1,intervalr=c(0.03,10)),spline.control = list(deg = 3,nknots = 10),lamv = seq(1,50,length.out = 20))

out.splogit.02<- link.compare.n1(model = 'splogit',ns = ns0,muv=-0.5,nrep = nrep0,model.args = list(r=0.2),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)),lamv = seq(1,50,length.out = 20))


out.splogit.15<- link.compare.n1(model = 'splogit',s0=0,ns = ns0,muv=-0.5,nrep = nrep0,model.args = list(r=1.5),init.args = list(init = c(0,0.2),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10)),bound=3,lamv = seq(5,50,length.out = 20))

rmse.out <-  cbind(out.logit$rmse.mat,out.probit$rmse.mat,out.robit3$rmse.mat,out.robit1$rmse.mat,out.robit2$rmse.mat, out.gev1$rmse.mat,out.gev2$rmse.mat,out.gev3$rmse.mat,out.gev4$rmse.mat,out.splogit.02$rmse.mat,out.splogit.15$rmse.mat)

prmse.out <-  cbind(out.logit$prmse.mat,out.probit$prmse.mat,out.robit3$prmse.mat,out.robit1$prmse.mat,out.robit2$prmse.mat, out.gev1$prmse.mat,out.gev2$prmse.mat,out.gev3$prmse.mat,out.gev4$prmse.mat,out.splogit.02$prmse.mat,out.splogit.15$prmse.mat)


col.name <-  c('logit','probit','robit','gev','splogit','pspline')
row.name <- c('logit','probit','robit(0.6)','robit(1)','robit(2)','gev(1)','gev(0.5)','gev(-0.5)',"gev(-1)",'splogit(0.2)','splogit(1.5)')

library(gridExtra)
save(rmse.out,prmse.out,file='output/output100_nonlinear.RData')
source('link_fun_code/tab.fig.fun.R')
res <- tab.fig.fun(rmse.out,col.name = col.name,row.name = row.name,remove = FALSE)
resp <- tab.fig.fun(prmse.out,col.name=col.name,row.name=row.name,remove=FALSE)
gg1 <- res$gp + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=11))
gg2 <- resp$gp + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=11))

pdf('figures/plot100n_rel.pdf',width = 12,height = 6)
grid.arrange(gg1,gg2,ncol=2)
dev.off()


######## nonlinear like linear -0.2(x-3)^2-0.15x+4 ######
source('link_fun_code/link.compare.n2.R')
ns0 <- 500
nrep0 <- 100

out.logit <- link.compare.n2(model = 'logit',ns = ns0,s0=0,nrep = nrep0,muv =-0.5,init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound = 3,lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.probit<- link.compare.n2(model = 'probit',ns = ns0,nrep = nrep0,muv =-0.5,init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound=3,lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.robit1<- link.compare.n2(model = 'robit',ns = ns0,nrep = nrep0,muv =-0.5,model.args = list(nu=1),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound=3,lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.robit2<- link.compare.n2(model = 'robit',ns = ns0,nrep = nrep0,muv =-0.5,model.args = list(nu=2),init.args = list(init = c(0,0.1),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound = 3,lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.robit3<- link.compare.n2(model = 'robit',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(nu=0.6),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound=3,lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.gev1 <- link.compare.n2(model = 'gev',ns = ns0,s0=0,nrep = nrep0,muv=-0.5,model.args = list(xi=1,locv=-1.5),init.args = list(init = c(0,0),xi0=0.6,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.gev2 <- link.compare.n2(model = 'gev',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(xi=0.5,locv=-1),init.args = list(init = c(0,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.gev3 <- link.compare.n2(model = 'gev',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(xi=-0.5,locv=-0.5),init.args = list(init = c(0,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.gev4 <- link.compare.n2(model = 'gev',ns = ns0,s0 = 0,muv=-0.5,nrep = nrep0,model.args = list(xi=-1,locv=0.5),init.args = list(init = c(0,0.2),xi0=-0.5,nu0=1,r0=1,intervalr=c(0.03,10),interval.nu = c(0.2,10)),spline.control = list(deg = 3,nknots = 10),lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.splogit.06<- link.compare.n2(model = 'splogit',s0=0,ns = ns0,muv=-0.5,nrep = nrep0,model.args = list(r=0.6),init.args = list(init = c(0.1,0.2),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'))


out.splogit.15<- link.compare.n2(model = 'splogit',s0=0,ns = ns0,muv=-0.5,nrep = nrep0,model.args = list(r=1.5),init.args = list(init = c(0.1,0.2),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound=3,lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'))


save(out.logit,out.probit,out.robit1,out.robit2,out.robit3,out.gev1,out.gev2,out.gev3,out.gev4,out.splogit.06,out.splogit.15,file='output/output500_nonlinear_case2.RData')


rmse.out <-  cbind(out.logit$rmse.mat,out.probit$rmse.mat,out.robit3$rmse.mat,out.robit1$rmse.mat,out.robit2$rmse.mat, out.gev1$rmse.mat,out.gev2$rmse.mat,out.gev3$rmse.mat,out.gev4$rmse.mat,out.splogit.06$rmse.mat,out.splogit.15$rmse.mat)

prmse.out <-  cbind(out.logit$prmse.mat,out.probit$prmse.mat,out.robit3$prmse.mat,out.robit1$prmse.mat,out.robit2$prmse.mat, out.gev1$prmse.mat,out.gev2$prmse.mat,out.gev3$prmse.mat,out.gev4$prmse.mat,out.splogit.06$prmse.mat,out.splogit.15$prmse.mat)


col.name <-  c('logit','probit','robit','gev','splogit','pspline')
row.name <- c('logit','probit','robit(0.6)','robit(1)','robit(2)','gev(1)','gev(0.5)','gev(-0.5)',"gev(-1)",'splogit(0.6)','splogit(1.5)')

library(gridExtra)
save(rmse.out,prmse.out,file='output/output500_nonlinear2.RData')
source('link_fun_code/tab.fig.fun.R')
res <- tab.fig.fun(rmse.out,col.name = col.name,row.name = row.name,remove = FALSE)
resp <- tab.fig.fun(prmse.out,col.name=col.name,row.name=row.name,remove=FALSE)
gg1 <- res$gp + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=11))
gg2 <- resp$gp + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=11))

pdf('figures/plot100n_rel.pdf',width = 12,height = 6)
grid.arrange(gg1,gg2,ncol=2)
dev.off()



###### case 3 ########
source('link_fun_code/link.compare.n2.R')
ns0 <- 500
nrep0 <- 100

out.logit <- link.compare.n2(model = 'logit',ns = ns0,s0=0,nrep = nrep0,muv =-0.5,init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.3,10)),bound = 3,lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)

out.probit<- link.compare.n2(model = 'probit',ns = ns0,nrep = nrep0,muv =-0.5,init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.4,10)),bound=3,lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)

out.robit1<- link.compare.n2(model = 'robit',ns = ns0,nrep = nrep0,muv =-0.5,model.args = list(nu=1),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound=3,lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)

out.robit2<- link.compare.n2(model = 'robit',ns = ns0,nrep = nrep0,muv =-0.5,model.args = list(nu=2),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.3,10)),bound = 3,lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)

out.robit3<- link.compare.n2(model = 'robit',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(nu=0.6),init.args = list(init = c(0,0.1),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound=3,lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)

out.gev1 <- link.compare.n2(model = 'gev',ns = ns0,s0=0,nrep = nrep0,muv=-0.5,model.args = list(xi=1,locv=-2),init.args = list(init = c(0,0),xi0=0.6,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)

out.gev2 <- link.compare.n2(model = 'gev',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(xi=0.5,locv=-1),init.args = list(init = c(0,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)

out.gev3 <- link.compare.n2(model = 'gev',ns = ns0,s0=11,nrep = nrep0,muv=-0.5,model.args = list(xi=-0.5,locv=1),init.args = list(init = c(0,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)

out.gev4 <- link.compare.n2(model = 'gev',ns = ns0,s0 = 0,muv=-0.5,nrep = nrep0,model.args = list(xi=-1,locv=1.5),init.args = list(init = c(0,0),xi0=-0.5,nu0=1,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),spline.control = list(deg = 3,nknots = 10),lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)


out.splogit.06<- link.compare.n2(model = 'splogit',s0=0,ns = ns0,muv=-0.5,nrep = nrep0,model.args = list(r=0.6),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)


out.splogit.15<- link.compare.n2(model = 'splogit',s0=0,ns = ns0,muv=-0.5,nrep = nrep0,model.args = list(r=1.5),init.args = list(init = c(0,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.5,10)),bound=3,lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)


save(out.logit,out.probit,out.robit1,out.robit2,out.robit3,out.gev1,out.gev2,out.gev3,out.gev4,out.splogit.06,out.splogit.15,file='output/output500_nonlinear_case3.RData')



###################case 1 ######################
source('link_fun_code/link.compare.n2.R')
ns0 <- 100
nrep0 <- 100

out.logit <- link.compare.n2(model = 'logit',ns = ns0,s0=0,nrep = nrep0,muv =-0.5,init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.3,10)),bound = 3,lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=1)

out.probit<- link.compare.n2(model = 'probit',ns = ns0,nrep = nrep0,muv =-0.5,init.args = list(init = c(0,0.1),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.5,10)),bound=3,lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=1)

out.robit1<- link.compare.n2(model = 'robit',ns = ns0,nrep = nrep0,muv =-0.5,model.args = list(nu=1),init.args = list(init = c(0,0.1),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound=3,lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=1)

out.robit2<- link.compare.n2(model = 'robit',ns = ns0,nrep = nrep0,muv =-0.5,model.args = list(nu=2),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.3,10)),bound = 3,lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=1)

out.robit3<- link.compare.n2(model = 'robit',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(nu=0.6),init.args = list(init = c(0,0.1),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound=3,lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=1)

out.gev1 <- link.compare.n2(model = 'gev',ns = ns0,s0=0,nrep = nrep0,muv=-0.5,model.args = list(xi=1,locv=-2),init.args = list(init = c(0,0),xi0=0.6,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=1)

out.gev2 <- link.compare.n2(model = 'gev',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(xi=0.5,locv=-1),init.args = list(init = c(0,0.1),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.3,10)),lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=1)

out.gev3 <- link.compare.n2(model = 'gev',ns = ns0,s0=11,nrep = nrep0,muv=-0.5,model.args = list(xi=-0.5,locv=1),init.args = list(init = c(0,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.3,10)),lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=1)

out.gev4 <- link.compare.n2(model = 'gev',ns = ns0,s0 = 0,muv=-0.5,nrep = nrep0,model.args = list(xi=-1,locv=1.5),init.args = list(init = c(0,0),xi0=-0.5,nu0=1,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),spline.control = list(deg = 3,nknots = 10),lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=1)


out.splogit.06<- link.compare.n2(model = 'splogit',s0=0,ns = ns0,muv=-0.5,nrep = nrep0,model.args = list(r=0.6),init.args = list(init = c(0,0.1),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.3,10)),lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=1)


out.splogit.15<- link.compare.n2(model = 'splogit',s0=0,ns = ns0,muv=-0.5,nrep = nrep0,model.args = list(r=1.5),init.args = list(init = c(0,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.5,10)),bound=3,lamv=seq(1,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=1)


save(out.logit,out.probit,out.robit1,out.robit2,out.robit3,out.gev1,out.gev2,out.gev3,out.gev4,out.splogit.06,out.splogit.15,file='output/output100_nonlinear_case1.RData')



######## algorithm in the paper #####
source('link_fun_code/link.compare.b4.R')
ns0 <- 100
nrep0 <- 100

out.logit <- link.compare.b4(model = 'logit',ns = ns0,nrep = nrep0,s0=0,
                             muv = -0.5,model.args = list(beta0=c(0,1,1)),
                             init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.1,10)),lamv=seq(2,200,length.out = 30), spline.control = list(deg = 3,nknots = 10,dd=1),weights.arg=c('equal','both','left','right'))

out.probit<- link.compare.b4(model = 'probit',ns = ns0,nrep = nrep0,muv = -0.5,model.args = list(beta0=c(0,1,1)),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=seq(2,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.robit1<- link.compare.b4(model = 'robit',ns = ns0,nrep = nrep0,muv = -0.5,model.args = list(beta0=c(0,1,1),nu=1),init.args = list(init = c(0,0,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=seq(2,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.robit2<- link.compare.b4(model = 'robit',ns = ns0,nrep = nrep0,muv = -0.5,model.args = list(beta0=c(0,1,1),nu=2),init.args = list(init = c(0,0,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound = 3,lamv=seq(2,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.robit3<- link.compare.b4(model = 'robit',ns = ns0,s0=0,nrep = nrep0,muv = -0.5,model.args = list(beta0=c(0,1,1),nu=0.6),init.args = list(init = c(0,0,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=seq(2,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.gev1 <- link.compare.b4(model = 'gev',ns = ns0,nrep =nrep0,muv = -0.5,s0=17,iter = 1000, model.args = list(beta0=c(0,1,1),xi=1,locv=-1.5),init.args = list(init = c(0,0,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound=3,spline.control = list(deg = 3,nknots = 10,dd=1),lamv=seq(2,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.gev2 <- link.compare.b4(model = 'gev',ns = ns0,nrep = nrep0,muv = -0.5,s0=0, model.args = list(beta0=c(0,1,1),xi=0.5,locv=-1),init.args = list(init = c(0,0,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),spline.control = list(deg = 3,nknots = 10,dd=1),lamv=seq(2,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.gev3 <- link.compare.b4(model = 'gev',ns = ns0,nrep = nrep0,muv = -0.5,model.args = list(beta0=c(0,1,1),xi=-0.5,locv=0),init.args = list(init = c(0,0.1,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),spline.control = list(deg = 3,nknots = 10,dd=1),lamv=seq(2,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.gev4 <- link.compare.b4(model = 'gev',ns = ns0,nrep = nrep0,muv = -0.5,s0=0,model.args = list(beta0=c(0,1,1),xi=-1,locv=1.2),init.args = list(init = c(0,0.1,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),spline.control = list(deg = 3,nknots = 10,dd=1),lamv=seq(2,200,length.out = 30),weights.arg=c('equal','both','left','right'))


out.splogit.02<- link.compare.b4(model = 'splogit',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(beta0=c(0,1,1),r=0.2),init.args = list(init = c(0,0,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.3,10)),lamv=seq(2,200,length.out = 30),weights.arg=c('equal','both','left','right'))


out.splogit.5<- link.compare.b4(model = 'splogit',s0=0,ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(beta0=c(0,1,1),r=5),init.args = list(init = c(0,0,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound = 3,lamv=seq(2,200,length.out = 30),weights.arg=c('equal','both','left','right'))


save(out.logit,out.probit,out.robit1,out.robit2,out.robit3,out.gev2,out.gev3,out.gev4,out.splogit.02,out.splogit.5,file='output/output100_binary_paper.RData')
