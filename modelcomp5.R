###  model compare new algorithm#### 
###### comparison of two covariates ####
source('link_fun_code/link.compare.b5.R')
ns0 <- 500
nrep0 <- 100

out.logit <- link.compare.b5(model = 'logit',ns = ns0,nrep = nrep0,s0=0,
                              muv = -0.5,model.args = list(beta0=c(0,1,1)),
                              init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.1,10)),lamv=10^(seq(-5,10,length.out = 50)), spline.control = list(deg = 3,nknots = 11,qv=1),weights.arg=c('equal','both','left','right'),iter=200)

out.probit<- link.compare.b5(model = 'probit',ns = ns0,nrep = nrep0,s0=0,
                             muv = -0.5,model.args = list(beta0=c(0,1,1),interval.nu=c(0.1,10)),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=10^(seq(-5,10,length.out = 50)), spline.control = list(deg = 3,nknots = 11,qv=1),weights.arg=c('equal','both','left','right'),iter=200)

out.robit1<- link.compare.b5(model = 'robit',ns = ns0,nrep = nrep0,muv = -0.5,model.args = list(beta0=c(0,1,1),nu=1),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=10^(seq(-5,10,length.out = 50)), spline.control = list(deg = 3,nknots = 11,qv=1),weights.arg=c('equal','both','left','right'),iter=200)

out.robit2<- link.compare.b5(model = 'robit',ns = ns0,nrep = nrep0,muv = -0.5,model.args = list(beta0=c(0,1,1),nu=2),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=10^(seq(-5,10,length.out = 50)), spline.control = list(deg = 3,nknots = 11,qv=1),weights.arg=c('equal','both','left','right'),iter=200)

out.robit3<- link.compare.b5(model = 'robit',ns = ns0,s0=0,nrep = nrep0,muv = -0.5,model.args = list(beta0=c(0,1,1),nu=0.6),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=10^(seq(-5,10,length.out = 50)), spline.control = list(deg = 3,nknots = 11,qv=1),weights.arg=c('equal','both','left','right'),iter=200)

out.gev1 <- link.compare.b5(model = 'gev',ns = ns0,nrep =nrep0,muv = -0.5,s0=0, model.args = list(beta0=c(0,1,1),xi=1,locv=-1.5),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=10^(seq(-5,10,length.out = 50)), spline.control = list(deg = 3,nknots = 11,qv=1),weights.arg=c('equal','both','left','right'),iter=200)

out.gev2 <- link.compare.b5(model = 'gev',ns = ns0,nrep = nrep0,muv = -0.5,s0=0, model.args = list(beta0=c(0,1,1),xi=0.5,locv=-1),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=10^(seq(-5,10,length.out = 50)), spline.control = list(deg = 3,nknots = 11,qv=1),weights.arg=c('equal','both','left','right'),iter=200)

out.gev3 <- link.compare.b5(model = 'gev',ns = ns0,nrep = nrep0,s0=0,muv = -0.5,model.args = list(beta0=c(0,1,1),xi=-0.5,locv=0),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=10^(seq(-5,10,length.out = 100)), spline.control = list(deg = 3,nknots = 11,qv=1),weights.arg=c('equal','both','left','right'),iter=200)

out.gev4 <- link.compare.b5(model = 'gev',ns = ns0,nrep = nrep0,muv = -0.5,s0=0,model.args = list(beta0=c(0,1,1),xi=-1,locv=1.2),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=10^(seq(-5,10,length.out = 100)), spline.control = list(deg = 3,nknots = 11,qv=1),weights.arg=c('equal','both','left','right'),iter=200)


out.splogit.02<- link.compare.b5(model = 'splogit',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(beta0=c(0,1,1),r=0.2),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.3,10)),lamv=10^(seq(-5,10,length.out = 50)), spline.control = list(deg = 3,nknots = 11,qv=1),weights.arg=c('equal','both','left','right'),iter=200)


out.splogit.5<- link.compare.b5(model = 'splogit',s0=0,ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(beta0=c(0,1,1),r=5),init.args = list(init = c(0,0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=10^(seq(-5,10,length.out = 50)), spline.control = list(deg = 3,nknots = 11,qv=1),weights.arg=c('equal','both','left','right'),iter=200)


save(out.logit,out.probit,out.robit3,out.robit1,out.robit2,out.gev1,out.gev2,out.gev3,out.gev4,out.splogit.02,out.splogit.5,file='output/new/output500_binary.RData')


######## nonlinear like linear -0.2(x-3)^2+4 case 2######
source('link_fun_code/link.compare.n5.R')
ns0 <- 100
nrep0 <- 100

out.logit <- link.compare.n5(model = 'logit',ns = ns0,s0=0,nrep = nrep0,muv =-0.5,init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.1,10)),bound = 3,lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'))

out.probit<- link.compare.n5(model = 'probit',ns = ns0,s0=0,nrep = nrep0,muv =-0.5,init.args = list(init = c(0,0.1),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound=3,lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'))

out.robit1<- link.compare.n5(model = 'robit',ns = ns0,nrep = nrep0,muv =-0.5,model.args = list(nu=1),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound=3,lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'))

out.robit2<- link.compare.n5(model = 'robit',ns = ns0,nrep = nrep0,muv =-0.5,model.args = list(nu=2),init.args = list(init = c(0,0.1),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound = 3,lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'))

out.robit3<- link.compare.n5(model = 'robit',ns = ns0,s0=0,nrep = nrep0,muv=-0.5,model.args = list(nu=0.6),init.args = list(init = c(0.1,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound=3,lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'))

out.gev1 <- link.compare.n5(model = 'gev',ns = ns0,s0=0,nrep = nrep0,muv=-0.5,model.args = list(xi=1,locv=-1.5),init.args = list(init = c(0,0),xi0=0.6,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.1,10)),lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'))

out.gev2 <- link.compare.n5(model = 'gev',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(xi=0.5,locv=-1),init.args = list(init = c(0,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.1,10)),lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'))

out.gev3 <- link.compare.n5(model = 'gev',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(xi=-0.5,locv=-0.5),init.args = list(init = c(0,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'))

out.gev4 <- link.compare.n5(model = 'gev',ns = ns0,s0 = 0,muv=-0.5,nrep = nrep0,model.args = list(xi=-1,locv=0.5),init.args = list(init = c(0,0.2),xi0=-0.5,nu0=1,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),spline.control = list(deg = 3,nknots = 11),lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'))

out.splogit.06<- link.compare.n5(model = 'splogit',s0=0,ns = ns0,muv=-0.5,nrep = nrep0,model.args = list(r=0.6),init.args = list(init = c(0.1,0.2),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.1,10)),lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'))


out.splogit.15<- link.compare.n5(model = 'splogit',s0=0,ns = ns0,muv=-0.5,nrep = nrep0,model.args = list(r=1.5),init.args = list(init = c(0.1,0.2),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound=3,lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'))


save(out.logit,out.probit,out.robit3,out.robit1,out.robit2,out.gev1,out.gev2,out.gev3,out.gev4,out.splogit.06,out.splogit.15,file='output/new/output100_nonlinear_case2.RData')



###### case 3 ########
source('link_fun_code/link.compare.n5.R')
ns0 <- 500
nrep0 <- 100

out.logit <- link.compare.n5(model = 'logit',ns = ns0,s0=0,nrep = nrep0,muv =-0.5,init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound = 3,lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'),case=3)

out.probit<- link.compare.n5(model = 'probit',ns = ns0,s0=0,nrep = nrep0,muv =-0.5,init.args = list(init = c(0,0.1),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound=3,lamv=10^(seq(-5,10,length.out = 50)),spline.control=list(deg=3,nknots=11),weights.arg=c('equal','both','left','right'),case=3)

out.robit1<- link.compare.n5(model = 'robit',ns = ns0,nrep = nrep0,muv =-0.5,model.args = list(nu=1),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound=3,lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'),case=3)

out.robit2<- link.compare.n5(model = 'robit',ns = ns0,nrep = nrep0,muv =-0.5,model.args = list(nu=2),init.args = list(init = c(0,0.1),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.25,10)),bound = 3,lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'),case=3)

out.robit3<- link.compare.n5(model = 'robit',ns = ns0,s0=0,nrep = nrep0,muv=-0.5,model.args = list(nu=0.6),init.args = list(init = c(0.1,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound=3,lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'),case=3)

out.gev1 <- link.compare.n5(model = 'gev',ns = ns0,s0=0,nrep = nrep0,muv=-0.5,model.args = list(xi=1,locv=-2),init.args = list(init = c(0,0),xi0=0.6,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'),case=3)

out.gev2 <- link.compare.n5(model = 'gev',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(xi=0.5,locv=-1),init.args = list(init = c(0,0),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'),case=3)

out.gev3 <- link.compare.n5(model = 'gev',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(xi=-0.5,locv=1),init.args = list(init = c(0,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'),case=3)

out.gev4 <- link.compare.n5(model = 'gev',ns = ns0,s0 = 0,muv=-0.5,nrep = nrep0,model.args = list(xi=-1,locv=1.5),init.args = list(init = c(0,0.2),xi0=-0.5,nu0=1,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),spline.control = list(deg = 3,nknots = 11),lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'),case=3)

out.splogit.06<- link.compare.n5(model = 'splogit',s0=0,ns = ns0,muv=-0.5,nrep = nrep0,model.args = list(r=0.6),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'),case=3)


out.splogit.15<- link.compare.n5(model = 'splogit',s0=0,ns = ns0,muv=-0.5,nrep = nrep0,model.args = list(r=1.5),init.args = list(init = c(0,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.5,10)),bound=3,lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'),case=3)


save(out.logit,out.probit,out.robit3,out.robit1,out.robit2,out.gev1,out.gev2,out.gev3,out.gev4,out.splogit.06,out.splogit.15,file='output/new/output100_nonlinear_case3.RData')


########### case 1 #################
source('link_fun_code/link.compare.n5.R')
ns0 <- 500 
nrep0 <- 100

out.logit <- link.compare.n5(model = 'logit',ns = ns0,s0=0,nrep = nrep0,muv =-0.5,init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.3,10)),bound = 3,lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'),case=1)

out.probit<- link.compare.n5(model = 'probit',ns = ns0,nrep = nrep0,muv =-0.5,init.args = list(init = c(0,0.1),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.5,10)),bound=3,lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'),case=1)

out.robit1<- link.compare.n5(model = 'robit',ns = ns0,nrep = nrep0,muv =-0.5,model.args = list(nu=1),init.args = list(init = c(0,0.1),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound=3,lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'),case=1)

out.robit2<- link.compare.n5(model = 'robit',ns = ns0,nrep = nrep0,muv =-0.5,model.args = list(nu=2),init.args = list(init = c(0,0),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.3,10)),bound = 3,lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'),case=1)

out.robit3<- link.compare.n5(model = 'robit',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(nu=0.6),init.args = list(init = c(0,0.1),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),bound=3,lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'),case=1)

out.gev1 <- link.compare.n5(model = 'gev',ns = ns0,s0=0,nrep = nrep0,muv=-0.5,model.args = list(xi=1,locv=-2),init.args = list(init = c(0,0),xi0=0.6,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'),case=1)

out.gev2 <- link.compare.n5(model = 'gev',ns = ns0,nrep = nrep0,muv=-0.5,model.args = list(xi=0.5,locv=-1),init.args = list(init = c(0,0.1),xi0=0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.3,10)),lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'),case=1)

out.gev3 <- link.compare.n5(model = 'gev',ns = ns0,s0=11,nrep = nrep0,muv=-0.5,model.args = list(xi=-0.5,locv=1),init.args = list(init = c(0,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.3,10)),lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'),case=1)

out.gev4 <- link.compare.n5(model = 'gev',ns = ns0,s0 = 0,muv=-0.5,nrep = nrep0,model.args = list(xi=-1,locv=1.5),init.args = list(init = c(0,0),xi0=-0.5,nu0=1,r0=1,intervalr=c(0.03,10),interval.nu=c(0.2,10)),spline.control = list(deg = 3,nknots = 10),lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'),case=1)


out.splogit.06<- link.compare.n5(model = 'splogit',s0=0,ns = ns0,muv=-0.5,nrep = nrep0,model.args = list(r=0.6),init.args = list(init = c(0,0.1),xi0=1,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.3,10)),lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'),case=1)


out.splogit.15<- link.compare.n5(model = 'splogit',s0=0,ns = ns0,muv=-0.5,nrep = nrep0,model.args = list(r=1.5),init.args = list(init = c(0,0),xi0=-0.5,nu0=2,r0=1,intervalr=c(0.03,10),interval.nu=c(0.5,10)),bound=3,lamv=10^(seq(-5,10,length.out = 50)),weights.arg=c('equal','both','left','right'),case=1)


save(out.logit,out.probit,out.robit1,out.robit2,out.robit3,out.gev1,out.gev2,out.gev3,out.gev4,out.splogit.06,out.splogit.15,file='output/new/output500_nonlinear_case1.RData')
