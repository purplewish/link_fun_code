##### based on output draw figures #####
library(reshape2)
library(ggplot2)
### remove means removing the gev results based on optim #####
tab.fig.fun <- function(mse.out,remove =TRUE)
{
  col.name <-  c('logit','probit','robit','gev0','gev','splogit','pspline')
  row.name <- c('logit','probit','robit(1)','robit(2)','gev(1)','gev(-1)','splogit(.01)','splogit(.05)','splogit(2)','splogit(5)')
  mse.mat <- sqrt(matrix(colMeans(mse.out),ncol=7,byrow=TRUE))
  colnames(mse.mat) <- col.name
  rownames(mse.mat) <- row.name
  
  if(remove == TRUE)
  {
    mse.mat <- mse.mat[,-4]
  }
  
  
  ##### sqrt_mse / min among different link functions #####
  mse.comp <- t(apply(mse.mat,1,function(x){x/min(x)}))
  mse.comp <- cbind(index = 1:nrow(mse.comp),mse.comp)
  mse.comp <- as.data.frame(mse.comp)
  mse.melt <- melt(mse.comp,value.name = 'rmse',id.vars = 'index',variable.name = 'model')
  
  gp <- ggplot(data = mse.melt,aes(x=index,y=rmse,color=model))+geom_line()+scale_x_discrete(labels=row.name) + theme_bw()+ylab('ratio') + xlab('assumption')
  
  res <- list(mse = mse.comp[,-1], gp=gp)
  return(res)
}

