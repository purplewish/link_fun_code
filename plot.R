load('output/output800.RData')
library(reshape2)
library(ggplot2)
mse.res <- sqrt(matrix(colMeans(mse.out),ncol=7,byrow=TRUE))
mse.res <- mse.res[,-4]
colnames(mse.res) <- c('logit','probit','robit','gev','splogit','pspline')
min.vector <- apply(mse.res,1,min)
mse.rel <- mse.res/min.vector
mse.res <- as.data.frame(cbind(index=1:10,mse.res))
mse.rel <- as.data.frame(cbind(index=1:10,mse.rel))
mse.melt <- melt(mse.res,id.vars='index',variable.name='model')
rel.melt <- melt(mse.rel,id.vars='index',variable.name='model')

pdf('document/figures/comparison/plot800.pdf',width = 10,height = 6)
qplot(x=index,y = value,color=model,data = mse.melt,geom = 'line')+scale_x_discrete(breaks=1:10, labels=c('logit','probit','robit(1)','robit(2)','gev(1)','gev(-1)','splogit(0.01)','splogit(0.05)','splogit(2)','splogit(5)'))
dev.off()

pdf('document/figures/comparison/plot500.rel.pdf',width = 10,height = 6)
qplot(x=index,y = value,color=model,data = rel.melt,geom = 'line')+scale_x_discrete(breaks=1:10, labels=c('logit','probit','robit(1)','robit(2)','gev(1)','gev(-1)','splogit(0.01)','splogit(0.05)','splogit(2)','splogit(5)'))
dev.off()


library(xtable)
print(xtable(mse.rel, caption='table', label='tab',digits=rep(4,8)), caption.placement="top",sanitize.colnames.function = identity,include.rownames=FALSE)
