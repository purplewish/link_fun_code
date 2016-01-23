################ application ####################
source('link_fun_code/psplinelink5.R')
source('link_fun_code/robit.em.R')
source('link_fun_code/splogit.mle.R')
source('link_fun_code/gev.mle.R')
library(mgcv)
weco <- read.table('data/weco.dat.txt',header=TRUE)


####### training set and test set ######
set.seed(7)
nr <- nrow(weco)
index <- sample(1:nr,round(nr/4))
test.dat <- weco[index,]
train.dat <- weco[-index,]


y0 <- train.dat$kwit
xmat0 <- train.dat[,c('sex','dex')]
newdata <- test.dat[,c('sex','dex')]
lam.m<- pspline.gcv5(y0 = train.dat$kwit,xmat = train.dat[,c('sex','dex')],qv=1,catv = 'sex',monotone = TRUE,nknots = 11,beta0 = c(1,1),MaxIter = 200,lamv=10^seq(-5,5,length=100))

pspline.fit.m <- psplinelink5(y0 = train.dat$kwit,xmat = train.dat[,c('sex','dex')],qv=1,catv = 'sex',monotone = TRUE,nknots = 11,beta0 = c(1,1),lambda=lam.m,MaxIter = 1000)

lam.wm<- pspline.gcv5(y0 = train.dat$kwit,xmat = train.dat[,c('sex','dex')],qv=1,catv = 'sex',monotone = FALSE,nknots = 11,beta0 = c(1,1),MaxIter = 200,lamv=10^seq(-5,5,length=100))

pspline.fit.wm <- psplinelink5(y0 = train.dat$kwit,xmat = train.dat[,c('sex','dex')],qv=1,catv = 'sex',monotone = FALSE,nknots = 11,beta0 = c(1,1),lambda=lam.wm,MaxIter = 1000)

init.value <- c(0,0,0)

logit.fit <- glm(kwit~sex+dex,train.dat,family = binomial(link='logit'))
probit.fit <- glm(kwit~sex+dex,train.dat, family=binomial(link='probit'))
robit.fit <- robit.pxem(y0 = y0,x0 = as.matrix(xmat0),beta0 = c(0,0,0),nu0 = 2,tol = 1e-3,interval.nu = c(0.1,10)) 
gev.fit<- gev.mle.new(y0 = y0,x0 = as.matrix(xmat0),par0 = c(0,0,-0.1,1),maxeval = 50000)
splogit.fit <- splogit.mle(y0 = y0,x0 = as.matrix(xmat0),par0 = init.value,intervalr = c(0,5))
gam.fit<- gam(kwit~sex+s(dex),train.dat,family = binomial(link = 'logit'))

logit.pred <- predict(logit.fit,newdata = as.data.frame(newdata),type = 'response')
probit.pred <- predict(probit.fit,newdata = as.data.frame(newdata),type = 'response')
robit.pred <- predict.pxem(est.obj = robit.fit,newdata = as.matrix(newdata))
gev.pred <- predict.gev.new(est.obj = gev.fit,newdata = as.matrix(newdata))
splogit.pred <- predict.splogit(est.obj = splogit.fit,newdata = as.matrix(newdata))
pspline.pred.m <- predict.pspline5(est.obj = pspline.fit.m,newdata = newdata)
pspline.pred.wm <- predict.pspline5(est.obj = pspline.fit.wm,newdata = newdata)
gam.pred <- predict(gam.fit,newdata = as.data.frame(newdata),type = 'response')



test.error <- function(test.y,pred.y)
{
  mean(test.y!=pred.y)
}

error.logit <- test.error(test.dat$kwit,1*(logit.pred>=0.5))
error.probit <- test.error(test.dat$kwit,1*(probit.pred>=0.5))
error.robit <- test.error(test.dat$kwit,1*(robit.pred>=0.5))
error.gev <- test.error(test.dat$kwit,1*(gev.pred>=0.5))
error.splogit <- test.error(test.dat$kwit,1*(splogit.pred>=0.5))
error.pspline.m <- test.error(test.dat$kwit,1*(pspline.pred.m>=0.5))
error.pspline.wm <- test.error(test.dat$kwit,1*(pspline.pred.wm>=0.5))
error.gam <- test.error(test.dat$kwit,1*(gam.pred>=0.5))

c(error.logit,error.probit,error.robit,error.gev ,error.splogit,error.pspline.wm,error.pspline.m ,error.gam)

dfp <- cbind(xmat0,logit.fit$fitted.value,probit.fit$fitted.values,robit.fit$fitted.values,gev.fit$fitted.values,splogit.fit$fitted.values,pspline.fit.wm$fitted.values,pspline.fit.m$fitted.values,gam.fit$fitted.values)

colnames(dfp) <- c('sex','dex','logit','probit','robit','GEV','splogit','pspline_wm','pspline_m','gam')
library(reshape2)
df.melt <- melt(dfp,value.name = 'prob',id.vars = c('sex','dex'),variable.name = 'link')
df.melt$link <- as.factor(df.melt$link)
library(ggplot2)

 ggplot(subset(df.melt,link%in%c("logit","robit","pspline_m")),aes(x=dex,y=prob,linetype=link,color=link))+geom_line()+facet_wrap(~sex)+theme_bw()+theme(strip.background = element_rect(fill='white'))


nll <- function(fitted.values,y0)
{
  -sum(y0*log(fitted.values)+(1-y0)*log(1-fitted.values))
}



 