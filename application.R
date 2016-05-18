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
set.seed(2016122)
# cross validation ###
## k = 5 #
index <- sample(1:nr,replace = FALSE)
nk <- round(nr/5)
error_mat <- matrix(0,5,8)
for(k in 1:5)
{
  index1 <- (k-1)*nk+1
  index2 <- k*nk
  if(k==5)
  {
    index2 <- nr
  }
  test.dat <- weco[index1:index2,]
  train.dat <- weco[-(index1:index2),]
  
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
  robit.fit <- robit.pxem(y0 = y0,x0 = as.matrix(xmat0),beta0 = c(0,0,0),nu0 = 2,tol = 1e-6,interval.nu = c(0.1,10)) 
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
  
  error_mat[k,] <- c(error.logit,error.probit,error.robit,error.gev ,error.splogit,error.pspline.wm,error.pspline.m ,error.gam)
}






#### all data 
lam.m0<- pspline.gcv5(y0 = weco$kwit,xmat = weco[,c('sex','dex')],qv=1,catv = 'sex',monotone = TRUE,nknots = 11,beta0 = c(1,1),MaxIter = 200,lamv=10^seq(-5,5,length=100))

pspline.fit.m0 <- psplinelink5(y0 = weco$kwit,xmat = weco[,c('sex','dex')],qv=1,catv = 'sex',monotone = TRUE,nknots = 11,beta0 = c(1,1),lambda=lam.m0,MaxIter = 1000)

lam.wm0<- pspline.gcv5(y0 = weco$kwit,xmat = weco[,c('sex','dex')],qv=1,catv = 'sex',monotone = FALSE,nknots = 11,beta0 = c(1,1),MaxIter = 200,lamv=10^seq(-5,5,length=100))

pspline.fit.wm0 <- psplinelink5(y0 = weco$kwit,xmat = weco[,c('sex','dex')],qv=1,catv = 'sex',monotone = FALSE,nknots = 11,beta0 = c(1,1),lambda=lam.wm0,MaxIter = 1000)

init.value <- c(0,0,0)

logit.fit0 <- glm(kwit~sex+dex,weco,family = binomial(link='logit'))
probit.fit0 <- glm(kwit~sex+dex,weco, family=binomial(link='probit'))
robit.fit0 <- robit.pxem(y0 = weco$kwit,x0 = as.matrix(weco[,c("sex","dex")]),beta0 = c(0,1,0),nu0 = 0.5,tol = 1e-6,interval.nu = c(0.1,10)) 
gev.fit0<- gev.mle.new(y0 = weco$kwit,x0 = as.matrix(weco[,c("sex","dex")]),par0 = c(0,0,-0.1,1),maxeval = 50000)
splogit.fit0 <- splogit.mle(y0 = weco$kwit,x0 = as.matrix(weco[,c("sex","dex")]),par0 = init.value,intervalr = c(0,5))
gam.fit0<- gam(kwit~sex+s(dex),weco,family = binomial(link = 'logit'))



dfp <- cbind(weco[,c("sex","dex","kwit")],logit.fit0$fitted.value,probit.fit0$fitted.values,robit.fit0$fitted.values,gev.fit0$fitted.values,splogit.fit0$fitted.values,pspline.fit.wm0$fitted.values,pspline.fit.m0$fitted.values,gam.fit0$fitted.values)

colnames(dfp) <- c('sex','dex',"kwit",'logit','probit','robit','GEV','splogit','pspline_wm','pspline_m','gam')
library(reshape2)
df.melt <- melt(dfp,value.name = 'prob',id.vars = c('sex','dex','kwit'),variable.name = 'link')
df.melt$link <- as.factor(df.melt$link)
df.melt$sex <- as.factor(df.melt$sex)
levels(df.melt$sex) <- c("sex 0"," sex 1")
library(ggplot2)

pdf("document/figures/application.pdf",width = 8,height = 4)
 ggplot(subset(df.melt,link%in%c("logit","robit","GEV","pspline_m")),aes(x=dex,y=prob,linetype=link,color=link))+geom_line()+geom_point(aes(x=dex,y=kwit),colour="black",size=1)+facet_wrap(~sex)+theme_bw()+theme(strip.background = element_rect(fill='white'))
dev.off()




 