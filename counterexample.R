

source('link_fun_code/psplinelink3.R')
source('link_fun_code/psplinelink4.R')
set.seed(18)
ns <- 100
betav <- c(0,1,1)
muv <- -0.5
sdv <- 1
bound <- 3
lowp <- muv - bound*sdv
upp <- muv + bound*sdv
nknots <- 10


x1 <- rnorm(2*ns,muv,sd = sdv)
x1 <- (x1[x1 >= lowp & x1 <= upp])[1:ns]
x2 <- rbinom(ns,size=1,prob=0.5)
eta0 <- cbind(1,x1,x2)%*%betav
prob0 <- 1-pgev(-eta0,loc = -1.5,scale = 1,shape = 1)
y0 <- rbinom(ns,size = 1,prob = prob0)
xmat0 <- cbind(x1,x2)

lamv <-seq(1,200,length.out=30)
beta00 <- coef(glm(y0~xmat0,family=binomial(link = 'logit')))
lam<- pspline.gcv4(y0 = y0,xmat = xmat0,monotone = TRUE,nknots = 10,beta0 = beta00 ,MaxIter = 1000,lamv = lamv)

pspline.fit <- psplinelink4(y0 = y0,xmat = xmat0,monotone = TRUE,nknots = nknots,beta0 = beta00,lambda=lam,MaxIter = iter)


lam3<- pspline.gcv3(y0 = y0,xmat = xmat0,qv=1,catv = 'x2',monotone = TRUE,nknots = nknots,beta0 = c(1,1),MaxIter = 1000,lamv = lamv,dd=1)

pspline.fit3 <- psplinelink3(y0 = y0,xmat = xmat0,qv=1,catv = 'x2',monotone = TRUE,nknots = nknots,beta0 = c(1,1),lambda=lam,MaxIter = 1000,dd=1)