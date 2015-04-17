##### simualte more than one covariates ####
set.seed(11)
ns <- 200
beta0 <- c(0,1,-0.5)
x1 <- runif(ns,min = -1.5,max = 3)
x2 <- runif(ns,min = 0,max = 3)
yita0 <- cbind(1,x1,x2)%*%beta0
prob0 <- exp(yita0)/(1+exp(yita0))
y0 <- rbinom(ns,size=1,prob=prob0)
summary(prob0)
res.glm <- glm(y0~x1+x2,family = 'binomial')
plot(res.glm$fitted.values,prob0)
abline(0,1)

sum((res.glm$fitted.values - prob0)^2)


#### gev link ####

set.seed(4)
xi <- 1
ns <- 500
beta0 <- c(0,1,-0.5)
x1 <- runif(ns,min = -8,max = 1)
x2 <- runif(ns,min = 0,max =2)
yita0 <- cbind(1,x1,x2)%*%beta0
prob0 <- 1-pgev(-yita0,loc = 0,scale = 1,shape = xi)
y0 <- rbinom(ns,size = 1,prob = prob0)
summary(prob0)
table(y0)

#### nolinear ####
set.seed(4)
ns <- 500
x1 <- runif(ns,min = -4,max = -1)
x2 <- runif(ns,min = 2 ,max = 4)
yita0 <- (x1+x2)^2+ x1
prob0 <- exp(yita0)/(1+exp(yita0))
y0 <- rbinom(ns,size = 1,prob = prob0)
summary(prob0)
table(y0)

##### categorical data ####

set.seed(4)
ns <- 500
x1 <- runif(ns,min = -3,max = 2)
x2 <- rbinom(ns,size=1,prob=0.5)
yita0 <- x1+x2
prob0 <- exp(yita0)/(1+exp(yita0))
y0 <- rbinom(ns,size = 1,prob = prob0)
summary(prob0)
table(y0)



source('link_fun_code/gev.mle.R')
source('link_fun_code/splogit.mle.R')
source('link_fun_code/psplinelink1.R')
library(mgcv)
res.gev <- gev.mle.new(x0 = cbind(x1,x2),y0 = y0,par0 = c(0,0.1,-0.5,-0.5),maxeval = 50000)

res.splogit <- splogit.mle(y0 = y0,x0 = cbind(x1,x2),par0 = c(0,0,0),intervalr = c(0.03,10))

res.sp <- psplinelink1(y0 = y0,xmat = cbind(x1,x2),qv=0.999,monotone = TRUE,nknots = 20,beta0 = c(1,-1),lambda=100)
res.add <- gam(y0~s(x1)+x2,family = binomial(),data = as.data.frame(cbind(y0,x1,x2)))
res.glm <- glm(y0~x1+x2,family = 'binomial')


plot(prob0,res.gev$fitted.values)
plot(prob0,res.add$fitted.values,col='red')
points(prob0,res.sp$fitted.values,col='blue')
points(prob0,res.glm$fitted.values,col='green')
abline(0,1)

sqrt(mean((res.glm$fitted.values - prob0)^2))
sqrt(mean((res.sp$fitted.values - prob0)^2))
sqrt(mean((res.splogit$fitted.values - prob0)^2))
sqrt(mean((res.add$fitted.values - prob0)^2))
sqrt(mean((res.gev$fitted.values - prob0)^2))



######might usefull ######
### standardized x and determine the value of a###
xmat <- cbind(x1,x2)
xmats <- scale(xmat)

atu <- quantile(sqrt(rowSums(xmats^2)),0.999)

### initial value of beta ####
glm.fit.logit <- glm(y0~0+x1+x2,data = as.data.frame(cbind(y0,xmats)),family = binomial())
beta00 <- coef(glm.fit.logit)
beta00 <- beta00/sqrt(sum(beta00^2))
beta00 <- sqrt(0.5)*sign(beta00)
eta.old <- xmats%*%beta00
q.old <- (eta.old/atu+1)/2
d.value <- 2
Ut <- pgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )

nknots <- 5
bs.nc <- nknots+deg-1
delta0 <- rep(0,bs.nc)
knots <- seq(0,1,length.out = nknots)

beta.old <- beta00
delta.old <- delta0

