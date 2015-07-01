################ application ####################
source('link_fun_code/psplinelink3.R')
source('link_fun_code/psplinelink4.R')
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
lam<- pspline.gcv3(y0 = train.dat$kwit,xmat = train.dat[,c('sex','dex')],qv=1,catv = 'sex',monotone = TRUE,nknots = 12,beta0 = c(1,1),MaxIter = 1000,lamv = seq(1,50,length.out=20),dd=1)

pspline.fit <- psplinelink3(y0 = train.dat$kwit,xmat = train.dat[,c('sex','dex')],qv=1,catv = 'sex',monotone = TRUE,nknots = 12,beta0 = c(1,1),lambda=,MaxIter = 1000,dd=1)

lam1<- pspline.gcv3(y0 = train.dat$kwit,xmat = train.dat[,c('sex','dex')],qv=1,catv = 'sex',monotone = FALSE,nknots = 12,beta0 = c(1,1),MaxIter = 1000,lamv = seq(1,50,length.out=20),dd=1)

pspline.fit1 <- psplinelink3(y0 = train.dat$kwit,xmat = train.dat[,c('sex','dex')],qv=1,catv = 'sex',monotone = FALSE,nknots = 12,beta0 = c(1,1),lambda=,MaxIter = 1000,dd=1)






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
pspline.pred <- predict.pspline3(est.obj = pspline.fit,newdata = newdata)
pspline.pred1 <- predict.pspline3(est.obj = pspline.fit1,newdata = newdata)
gam.pred <- predict(gam.fit,newdata = as.data.frame(newdata),type = 'response')



test.error <- function(test.y,pred.y)
{
  mean(test.y!=pred.y)
}
# table(test.dat$kwit,1*(logit.pred>=0.5))
# table(test.dat$kwit,1*(probit.pred>=0.5))
# table(test.dat$kwit,1*(robit.pred>=0.5))
# table(test.dat$kwit,1*(gev.pred>=0.5))
# table(test.dat$kwit,1*(splogit.pred>=0.5))
# table(test.dat$kwit,1*(pspline.pred>=0.5))
# table(test.dat$kwit,1*(pspline.pred1>=0.5))
# table(test.dat$kwit,1*(gam.pred>=0.5))


error.logit <- test.error(test.dat$kwit,1*(logit.pred>=0.5))
error.probit <- test.error(test.dat$kwit,1*(probit.pred>=0.5))
error.robit <- test.error(test.dat$kwit,1*(robit.pred>=0.5))
error.gev <- test.error(test.dat$kwit,1*(gev.pred>=0.5))
error.splogit <- test.error(test.dat$kwit,1*(splogit.pred>=0.5))
error.pspline <- test.error(test.dat$kwit,1*(pspline.pred>=0.5))
error.pspline1 <- test.error(test.dat$kwit,1*(pspline.pred1>=0.5))
error.gam <- test.error(test.dat$kwit,1*(gam.pred>=0.5))

dfp <- cbind(xmat0,logit.fit$fitted.value,probit.fit$fitted.values,robit.fit$fitted.values,gev.fit$fitted.values,splogit.fit$fitted.values,pspline.fit$fitted.values,pspline.fit1$fitted.values,gam.fit$fitted.values)

colnames(dfp) <- c('sex','dex','logit','probit','robit','GEV','splogit','pspline_m','pspline_wom','gam')
library(reshape2)
df.melt <- melt(dfp,value.name = 'prob',id.vars = c('sex','dex'),variable.name = 'link')
df.melt$link <- as.factor(df.melt$link)
library(ggplot2)

# ggplot(df.melt,aes(x=dex,y=prob,linetype=link,color=link))+geom_line()+facet_wrap(~sex)+theme_bw()+theme(strip.background = element_rect(fill='white'))


nll <- function(fitted.values,y0)
{
  -sum(y0*log(fitted.values)+(1-y0)*log(1-fitted.values))
}

# qchisq(p=0.95,df = pspline.fit$traceH)
# qchisq(p=0.95,df = 8+pspline.fit$traceH-1)
# 2*nll(gev.fit$fitted.values,y0)-2*nll(pspline.fit$fitted.values,y0)


beta0 <- coef(glm(kwit~sex+dex,weco,family = binomial(link='logit')))
 lam4<- pspline.gcv4(y0 = weco$kwit,xmat0 = weco[,c('sex','dex')],monotone = TRUE,nknots = 11,beta0 = beta0,MaxIter = 1000,lamv = seq(1,50,length.out=20))
# 
 pspline.fit4 <- psplinelink4(y0 = weco$kwit,xmat0 = weco[,c('sex','dex')],monotone = TRUE,nknots = 11,beta0 =beta0 ,lambda=lam4,MaxIter = 1000)

 lam3<- pspline.gcv3(y0 = weco$kwit,xmat = weco[,c('sex','dex')],qv=1,catv = 'sex',monotone = TRUE,nknots = 12,beta0 = c(1,1),MaxIter = 1000,lamv = seq(1,50,length.out=20),dd=1)
 
 pspline.fit3 <- psplinelink3(y0 = weco$kwit,xmat = weco[,c('sex','dex')],qv=1,catv = 'sex',monotone = TRUE,nknots = 12,beta0 = c(1,1),lambda=lam3,MaxIter = 1000,dd=1)
 


 