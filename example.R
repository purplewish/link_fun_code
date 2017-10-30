
source('link_fun_code/psplinelink4.R')
source("link_fun_code/psplinelink5.R")
source("link_fun_code/psplinelink3_v2.R")
library(evd)
library(xtable)
set.seed(20)
ns <- 100
betav <- c(0,1,1)
muv <- -0.5
sdv <- 1
bound <- 3
lowp <- muv - bound*sdv
upp <- muv + bound*sdv
nknots <- 11 ## including noninterior knots


x1 <- rnorm(2*ns,muv,sd = sdv)
x1 <- (x1[x1 >= lowp & x1 <= upp])[1:ns]
x2 <- rbinom(ns,size=1,prob=0.5)
eta0 <- cbind(1,x1,x2)%*%betav
prob0 <- 1-pgev(-eta0,loc = -1.5,scale = 1,shape = 1)
y0 <- rbinom(ns,size = 1,prob = prob0)
xmat0 <- cbind(x1,x2)


lamv=10^(seq(0,10,length.out = 20))
beta00 <- coef(glm(y0~xmat0,family=binomial(link = 'logit')))
beta01 <- coef(glm(y0~0+xmat0,family=binomial(link = 'logit')))
beta01 <- beta01/sqrt(sum(beta01^2))

#lam4 <- pspline.gcv4(y0 = y0,xmat0 = xmat0,monotone = TRUE,nknots = 10,beta0=beta00,lamv = lamv)

fit41 <- psplinelink4(y0 = y0,xmat0 = xmat0,monotone = TRUE,nknots = 11,beta0 = beta00,lambda=20,MaxIter = 1000,kp=1e6)
diff(fit41$delta)

fit42 <- psplinelink4(y0 = y0,xmat0 = xmat0,monotone = TRUE,nknots = 11,beta0 = beta00,lambda=20,MaxIter = 1000,kp=1e7)
diff(fit42$delta)

psplinelink4(y0 = y0,xmat0 = xmat0,monotone = TRUE,nknots = 11,beta0 = beta00,lambda=10,MaxIter = 1000,kp=1e14)



fit43 <- psplinelink4(y0 = y0,xmat0 = xmat0,monotone = TRUE,nknots = 11,beta0 = beta00,lambda=100,MaxIter = 1000,kp=1e6)
diff(fit43$delta)


tab1 <- t(cbind(diff(fit41$delta),diff(fit42$delta)))
print(xtable(rbind(tab1[,1:6],tab1[,7:12]),digits=-2))


set.seed(1050)
a0 <- 4
x11 <- rnorm(a0*ns,muv,sd = sdv)
x11 <- (x1[x1 >= lowp & x1 <= upp])[1:ns]
x21 <- rbinom(a0*ns,size=1,prob=0.5)

eta01 <- cbind(1,x11,x21)%*%betav
prob01 <- 1-pgev(-eta01,loc = -1.5,scale = 1,shape = 1)
y01 <- rbinom(a0*ns,size = 1,prob = prob01)
xmat01 <- rbind(xmat0,cbind(x11,x21))
y01 <- c(y0,y01)


fit44 <- psplinelink4(y0 = y01,xmat0 = xmat01,monotone = TRUE,nknots = 11,beta0 = beta00,lambda=20,MaxIter = 1000,kp=1e6)
diff(fit44$delta)


### compare two algorithm with different ####

lamv=10^(seq(-5,12,length.out = 100))
coefmat3 <- coefmat5 <- matrix(0,100,2)
for(i in 1:100)
{
  set.seed(18+i)
  ns <- 100
  betav <- c(0,1,1)
  muv <- -0.5
  sdv <- 1
  bound <- 3
  lowp <- muv - bound*sdv
  upp <- muv + bound*sdv
  nknots <- 11
  
  
  x1 <- rnorm(2*ns,muv,sd = sdv)
  x1 <- (x1[x1 >= lowp & x1 <= upp])[1:ns]
  x2 <- rbinom(ns,size=1,prob=0.5)
  eta0 <- cbind(1,x1,x2)%*%betav
  prob0 <- 1-pgev(-eta0,loc = -1.5,scale = 1,shape = 1)
  y0 <- rbinom(ns,size = 1,prob = prob0)
  xmat0 <- cbind(x1,x2)
  
  beta01 <- coef(glm(y0~0+xmat0,family=binomial(link = 'logit')))
  beta01 <- beta01/sqrt(sum(beta01^2))
  
  lam5 <- pspline.gcv5(y0 = y0,xmat = xmat0,qv=1,catv = 'x2',monotone = TRUE,nknots = 11,beta0 = beta01,MaxIter = 1000,lamv = lamv)
  
  fit5 <- psplinelink5(y0 = y0,xmat = xmat0,qv=1, catv ="x2",monotone = TRUE,nknots=11,beta0 = beta01,lambda = lam5)
  coefmat5[i,] <- fit5$est
  
  
  tryCatch({coefmat3[i,] <- psplinelink3(y0 = y0,xmat = xmat0,qv=1, catv ="x2",monotone = TRUE,nknots =11, beta0 = beta01,lambda = lam5)$est},
           error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


index3 <- which(rowSums(coefmat3)==0)


summary(c(coefmat3[-index3,] - coefmat5[-index3,]))

save(coefmat3,coefmat5,file = "document/STAT_COM/coef35.RData")


set.seed(18)
ns <- 100
betav <- c(0,1,1)
muv <- -0.5
sdv <- 1
bound <- 3
lowp <- muv - bound*sdv
upp <- muv + bound*sdv
nknots <- 11 ## including noninterior knots


x1 <- rnorm(2*ns,muv,sd = sdv)
x1 <- (x1[x1 >= lowp & x1 <= upp])[1:ns]
x2 <- rbinom(ns,size=1,prob=0.5)
eta0 <- cbind(1,x1,x2)%*%betav
prob0 <- 1-pgev(-eta0,loc = -1.5,scale = 1,shape = 1)
y0 <- rbinom(ns,size = 1,prob = prob0)
xmat0 <- cbind(x1,x2)


lamv=10^(seq(0,10,length.out = 20))
beta00 <- coef(glm(y0~xmat0,family=binomial(link = 'logit')))
beta01 <- coef(glm(y0~0+xmat0,family=binomial(link = 'logit')))
beta01 <- beta01/sqrt(sum(beta01^2))


fit45 <- psplinelink4(y0 = y0,xmat0 = xmat0,monotone = TRUE,nknots = 11,beta0 = beta00,lambda=1,MaxIter = 1000,kp=1e6)
