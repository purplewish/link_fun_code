### 
source('link_fun_code/robit.em.R')
source('link_fun_code/splogit.mle.R')
source('link_fun_code/splinelink.R')
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


zz <- t(apply(out.gev1$aic[,-4],1,function(x){x - x[4]}))

zznew <- cbind(index= 1:100,zz[,-4])
zznew <- as.data.frame(zznew)

zz.melt <- melt(data = zznew,id.vars = 'index',value.name = 'aic',variable.name='model')
g1 <- ggplot(zz.melt,aes(x=index,y=aic,color=model))+geom_point()+theme_bw()+ggtitle('gev')+geom_hline(yintercept=0)

zz1 <- t(apply(out.logit$aic[,-4],1,function(x){x - x[1]}))

zznew1 <- cbind(index= 1:100,zz1[,-1])
zznew1 <- as.data.frame(zznew1)

zz.melt1 <- melt(data = zznew1,id.vars = 'index',value.name = 'aic',variable.name='model')
g2 <- ggplot(zz.melt1,aes(x=index,y=aic,color=model))+geom_point()+theme_bw()+ggtitle('logit')+geom_hline(yintercept=0)



g12 <- list()
g12[[1]] <- g1
g12[[2]] <- g2

pdf('figures/gev_logit.pdf')
print(g12)
dev.off()
