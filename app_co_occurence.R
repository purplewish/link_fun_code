

###### co-occurence #####
setwd("~/Research/link_function/")
protea.dat <- read.csv('data/Jiang2013_data/ProteaData.csv')
#hist(protea.dat$observed/protea.dat$n.sites,breaks = seq(0,1,length.out = 21))
source('link_fun_code/psplinelink5.R')
source('link_fun_code/robit.em.R')
source('link_fun_code/splogit.mle.R')
source('link_fun_code/gev.mle.R')

#### 10 fold cross validation #### 

set.seed(7)
nr <- nrow(protea.dat)

covariate.names <- c("dist","FireSurv","Height","FlowDif","Pollen","LenWid","lfarea","p.expected.bayes")

cat.names <- c('FireSurv','FlowDif','Pollen')
protea.dat$p.expected.bayes <- log(protea.dat$p.expected.bayes)-log(1-protea.dat$p.expected.bayes)

set.seed(77)
index0 <- sample.int(nr)
subnr <- round(nr/10)
logit.error <- probit.error <- robit.error <- splogit.error <- gev.error <- pspline.error.m <- pspline.error.wm <- gam.error <- rep(0,10)
for(j in 1:10)
{
  if(j < 10) 
    {
    index1 <- (j-1)*subnr+1
    index2 <- subnr*j
    }
  
  if(j==10)
  {
    index1 <- 9*subnr+1
    index2 <- nr
  }
  train.dat <- protea.dat[-(index1:index2),c(covariate.names,"observed","n.sites")]
  test.dat <- protea.dat[index1:index2,c(covariate.names,"observed","n.sites")]
  
  xmat0 <- train.dat[,covariate.names]
  newdata <- test.dat[,covariate.names]
  logit.fit <- glm(cbind(observed,n.sites-observed)~.,train.dat,family = binomial(link='logit'))
  probit.fit <- glm(cbind(observed,n.sites-observed)~.,train.dat, family=binomial(link='probit'))
  robit.fit <- robit.pxem.binom(y0 = train.dat$observed,size=train.dat$n.sites,x0 = as.matrix(xmat0),beta0 = rep(0,9),nu0 = 2,tol = 1e-4,interval.nu=c(0.1,10)) 
  # gev.fit<- gev.mle.new(y0 = train.dat$observed,x0 = as.matrix(xmat0),size=train.dat$n.sites,par0 = c(rep(0,9),-0.5),maxeval = 50000)
  splogit.fit <- splogit.mle(y0 = train.dat$observed,size=train.dat$n.sites,x0 = as.matrix(xmat0),par0 = rep(0,9),intervalr = c(0,10))
  gam.fit<- gam(cbind(observed,n.sites)~FireSurv+FlowDif+Pollen+s(dist)+s(Height)+s(lfarea)+s(LenWid)+s(p.expected.bayes),train.dat,family = binomial(link = 'logit'))
  
  
  lam.m<- pspline.gcv5(y0 = train.dat$observed,xmat = as.matrix(xmat0),size = train.dat$n.sites,qv=1,catv = cat.names,monotone = TRUE,nknots = 20,beta0 = rep(1,8),MaxIter = 500,lamv = 10^seq(-5,10,length=100))
  
  pspline.fit.m <- psplinelink5(y0 = train.dat$observed,xmat = as.matrix(xmat0),size = train.dat$n.sites,qv=1,catv = cat.names,monotone = TRUE,nknots = 20,beta0 = rep(1,8),MaxIter = 500,lambda = lam.m)

  lam.wm<- pspline.gcv5(y0 = train.dat$observed,xmat = as.matrix(xmat0),size = train.dat$n.sites,qv=1,catv = cat.names,monotone = FALSE,nknots = 20,beta0 = rep(1,8),MaxIter = 500,lamv = 10^seq(-5,10,length=100))
  
  pspline.fit.wm <- psplinelink5(y0 = train.dat$observed,xmat = as.matrix(xmat0),size = train.dat$n.sites,qv=1,catv = cat.names,monotone = FALSE,nknots = 20,beta0 = rep(1,8),MaxIter = 500,lambda = lam.wm)
  

  logit.pred <- predict(logit.fit,newdata = as.data.frame(newdata),type = 'response')
  probit.pred <- predict(probit.fit,newdata = as.data.frame(newdata),type = 'response')
  robit.pred <- predict.pxem(est.obj = robit.fit,newdata = as.matrix(newdata))
  #gev.pred <- predict.gev.new(est.obj = gev.fit,newdata = as.matrix(newdata))
  splogit.pred <- predict.splogit(est.obj = splogit.fit,newdata = as.matrix(newdata))
  pspline.pred.m <- predict.pspline5(est.obj = pspline.fit.m,newdata = as.matrix(newdata))
  pspline.pred.wm <- predict.pspline5(est.obj = pspline.fit.wm,newdata = as.matrix(newdata))
  gam.pred <- predict(gam.fit,newdata = as.data.frame(newdata),type = 'response')
  
  logit.error[j] <- sum((logit.pred*test.dat$n.sites - test.dat$observed)^2)
  robit.error[j] <- sum((robit.pred*test.dat$n.sites - test.dat$observed)^2)
  probit.error[j] <- sum((probit.pred*test.dat$n.sites - test.dat$observed)^2)
  splogit.error[j] <- sum((splogit.pred*test.dat$n.sites - test.dat$observed)^2)
  #gev.error[j] <- sum((gev.pred*test.dat$n.sites - test.dat$observed)^2)
  gam.error[j] <- sum((gam.pred*test.dat$n.sites - test.dat$observed)^2)
  pspline.error.m[j] <- sum((pspline.pred.m*test.dat$n.sites - test.dat$observed)^2)
  pspline.error.wm[j] <- sum((pspline.pred.wm*test.dat$n.sites - test.dat$observed)^2)
  
  print(j)
  gc()
}


error <- c(sum(logit.error)/nr,sum(probit.error)/nr,sum(robit.error)/nr,sum(splogit.error)/nr,sum(gam.error)/nr,sum(pspline.error.wm)/nr,sum(pspline.error.m)/nr)

save(error,file='output/protea_error.RData')





