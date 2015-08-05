###  model compare new algorithm#### 
###### comparison of two covariates ####
source('link_fun_code/link.compare.b5.R')
ns0 <- 100
nrep0 <- 100

out.logit <- link.compare.b5(model = 'logit',ns = ns0,nrep = nrep0,s0=0,
                              muv = -0.5,model.args = list(beta0=c(0,1,1)),
                              init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.1,10)),lamv=10^(seq(-5,10,length.out = 50)), spline.control = list(deg = 3,nknots = 11,qv=0.95),weights.arg=c('equal','both','left','right'),iter=200)

out.probit<- link.compare.b5(model = 'probit',ns = ns0,nrep = nrep0,muv = -0.5,model.args = list(beta0=c(0,1,1),interval.nu=c(0.1,10)),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.1,10)),lamv=10^(seq(-5,10,length.out = 50)), spline.control = list(deg = 3,nknots = 11,qv=0.95),weights.arg=c('equal','both','left','right'),iter=200)

out.robit1<- link.compare.b5(model = 'robit',ns = ns0,nrep = nrep0,muv = -0.5,model.args = list(beta0=c(0,1,1),nu=1),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.1,10)),lamv=10^(seq(-5,10,length.out = 50)), spline.control = list(deg = 3,nknots = 11,qv=0.95),weights.arg=c('equal','both','left','right'),iter=200)

out.robit2<- link.compare.b5(model = 'robit',ns = ns0,nrep = nrep0,muv = -0.5,model.args = list(beta0=c(0,1,1),nu=2),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.1,10)),lamv=10^(seq(-5,10,length.out = 50)), spline.control = list(deg = 3,nknots = 11,qv=0.95),weights.arg=c('equal','both','left','right'),iter=200)

out.robit3<- link.compare.b5(model = 'robit',ns = ns0,s0=0,nrep = nrep0,muv = -0.5,model.args = list(beta0=c(0,1,1),nu=0.6),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.1,10)),lamv=10^(seq(-5,10,length.out = 50)), spline.control = list(deg = 3,nknots = 11,qv=0.95),weights.arg=c('equal','both','left','right'),iter=200)

out.gev1 <- link.compare.b5(model = 'gev',ns = ns0,nrep =nrep0,muv = -0.5,s0=0, model.args = list(beta0=c(0,1,1),xi=1,locv=-1.5),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.1,10)),lamv=10^(seq(-5,10,length.out = 50)), spline.control = list(deg = 3,nknots = 11,qv=0.95),weights.arg=c('equal','both','left','right'),iter=200)

out.gev2 <- link.compare.b5(model = 'gev',ns = ns0,nrep = nrep0,muv = -0.5,s0=0, model.args = list(beta0=c(0,1,1),xi=0.5,locv=-1),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.1,10)),lamv=10^(seq(-5,10,length.out = 50)), spline.control = list(deg = 3,nknots = 11,qv=0.95),weights.arg=c('equal','both','left','right'),iter=200)

out.gev3 <- link.compare.b5(model = 'gev',ns = ns0,nrep = nrep0,muv = -0.5,model.args = list(beta0=c(0,1,1),xi=-0.5,locv=0),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.1,10)),lamv=10^(seq(-5,10,length.out = 50)), spline.control = list(deg = 3,nknots = 11,qv=0.95),weights.arg=c('equal','both','left','right'),iter=200)

out.gev4 <- link.compare.b5(model = 'gev',ns = ns0,nrep = nrep0,muv = -0.5,s0=0,model.args = list(beta0=c(0,1,1),xi=-1,locv=1.2),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.1,10)),lamv=10^(seq(-5,10,length.out = 50)), spline.control = list(deg = 3,nknots = 11,qv=0.95),weights.arg=c('equal','both','left','right'),iter=200)


out.splogit.02<- link.compare.b5(model = 'splogit',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(beta0=c(0,1,1),r=0.2),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.1,10)),lamv=10^(seq(-5,10,length.out = 50)), spline.control = list(deg = 3,nknots = 11,qv=0.95),weights.arg=c('equal','both','left','right'),iter=200)


out.splogit.5<- link.compare.b5(model = 'splogit',s0=0,ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(beta0=c(0,1,1),r=5),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.1,10)),lamv=10^(seq(-5,10,length.out = 50)), spline.control = list(deg = 3,nknots = 11,qv=0.95),weights.arg=c('equal','both','left','right'),iter=200)


prmse.out <-  cbind(out.logit$prmse.mat,out.probit$prmse.mat,out.robit3$prmse.mat,out.robit1$prmse.mat,out.robit2$prmse.mat, out.gev1$prmse.mat,out.gev2$prmse.mat,out.gev3$prmse.mat,out.gev4$prmse.mat,out.splogit.02$prmse.mat, out.splogit.5$prmse.mat)

wprmse.out <-  cbind(out.logit$wprmse.mat,out.probit$wprmse.mat,out.robit3$wprmse.mat,out.robit1$wprmse.mat,out.robit2$wprmse.mat, out.gev1$wprmse.mat,out.gev2$wprmse.mat,out.gev3$wprmse.mat,out.gev4$wprmse.mat,out.splogit.02$wprmse.mat, out.splogit.5$wprmse.mat)


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


######## nonlinear like linear -0.2(x-3)^2+4 ######
source('link_fun_code/link.compare.n2.R')
ns0 <- 500
nrep0 <- 100

out.logit <- link.compare.n2(model = 'logit',ns = ns0,s0=0,nrep = nrep0,muv =-0.5,init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)),bound = 3,lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.probit<- link.compare.n2(model = 'probit',ns = ns0,nrep = nrep0,muv =-0.5,init.args = list(init = c(0,0.1),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)),bound=3,lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.robit1<- link.compare.n2(model = 'robit',ns = ns0,nrep = nrep0,muv =-0.5,model.args = list(nu=1),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)),bound=3,lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.robit2<- link.compare.n2(model = 'robit',ns = ns0,nrep = nrep0,muv =-0.5,model.args = list(nu=2),init.args = list(init = c(0,0.1),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)),bound = 3,lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.robit3<- link.compare.n2(model = 'robit',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(nu=0.6),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)),bound=3,lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.gev1 <- link.compare.n2(model = 'gev',ns = ns0,s0=0,nrep = nrep0,muv=-0.5,model.args = list(xi=1,locv=-2),init.args = list(init = c(0,0),xi0=0.6,nu0=2,r0=1,intervalr=c(0.03,10)),lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.gev2 <- link.compare.n2(model = 'gev',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(xi=0.5,locv=-1),init.args = list(init = c(0,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10)),lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.gev3 <- link.compare.n2(model = 'gev',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(xi=-0.5,locv=0),init.args = list(init = c(0,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10)),lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.gev4 <- link.compare.n2(model = 'gev',ns = ns0,s0 = 0,muv=-0.5,nrep = nrep0,model.args = list(xi=-1,locv=0),init.args = list(init = c(0,0.2),xi0=-0.5,nu0=1,r0=1,intervalr=c(0.03,10)),spline.control = list(deg = 3,nknots = 10),lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'))

out.splogit.06<- link.compare.n2(model = 'splogit',s0=0,ns = ns0,muv=-0.5,nrep = nrep0,model.args = list(r=0.6),init.args = list(init = c(0.1,0.2),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)),lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'))


out.splogit.15<- link.compare.n2(model = 'splogit',s0=0,ns = ns0,muv=-0.5,nrep = nrep0,model.args = list(r=1.5),init.args = list(init = c(0.1,0.2),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10)),bound=3,lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'))


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

out.logit <- link.compare.n2(model = 'logit',ns = ns0,s0=0,nrep = nrep0,muv =-0.5,init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)),bound = 3,lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)

out.probit<- link.compare.n2(model = 'probit',ns = ns0,nrep = nrep0,muv =-0.5,init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)),bound=3,lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)

out.robit1<- link.compare.n2(model = 'robit',ns = ns0,nrep = nrep0,muv =-0.5,model.args = list(nu=1),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)),bound=3,lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)

out.robit2<- link.compare.n2(model = 'robit',ns = ns0,nrep = nrep0,muv =-0.5,model.args = list(nu=2),init.args = list(init = c(0,0.1),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)),bound = 3,lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)

out.robit3<- link.compare.n2(model = 'robit',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(nu=0.6),init.args = list(init = c(0,0.1),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)),bound=3,lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)

out.gev1 <- link.compare.n2(model = 'gev',ns = ns0,s0=0,nrep = nrep0,muv=-0.5,model.args = list(xi=1,locv=-2),init.args = list(init = c(0,0),xi0=0.6,nu0=2,r0=1,intervalr=c(0.03,10)),lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)

out.gev2 <- link.compare.n2(model = 'gev',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(xi=0.5,locv=-1),init.args = list(init = c(0,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10)),lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)

out.gev3 <- link.compare.n2(model = 'gev',ns = ns0,s0=11,nrep = nrep0,muv=-0.5,model.args = list(xi=-0.5,locv=1),init.args = list(init = c(0,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10)),lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)

out.gev4 <- link.compare.n2(model = 'gev',ns = ns0,s0 = 0,muv=-0.5,nrep = nrep0,model.args = list(xi=-1,locv=1.5),init.args = list(init = c(0,0),xi0=-0.5,nu0=1,r0=1,intervalr=c(0.03,10)),spline.control = list(deg = 3,nknots = 10),lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)


out.splogit.06<- link.compare.n2(model = 'splogit',s0=0,ns = ns0,muv=-0.5,nrep = nrep0,model.args = list(r=0.6),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10)),lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)


out.splogit.15<- link.compare.n2(model = 'splogit',s0=0,ns = ns0,muv=-0.5,nrep = nrep0,model.args = list(r=1.5),init.args = list(init = c(0,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10)),bound=3,lamv=seq(5,200,length.out = 30),weights.arg=c('equal','both','left','right'),case=3)


save(out.logit,out.probit,out.robit1,out.robit2,out.robit3,out.gev1,out.gev2,out.gev3,out.gev4,out.splogit.06,out.splogit.15,file='output/output200_nonlinear_case3.RData')
