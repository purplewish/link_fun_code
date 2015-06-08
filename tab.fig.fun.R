##### based on output draw figures #####
library(reshape2)
library(ggplot2)
### remove means removing the gev results based on optim #####
#col.name <-  c('logit','probit','robit','gev','splogit','pspline')
#row.name <- c('logit','probit','robit(1)','robit(2)','gev(1)','gev(-1)','splogit(.01)','splogit(.05)','splogit(2)','splogit(5)')

tab.fig.fun <- function(rmse.out,remove =TRUE,col.name, row.name)
{ 
  rmse.mat <- matrix(colMeans(rmse.out),ncol=length(col.name),byrow=TRUE)
  colnames(rmse.mat) <- col.name
  rownames(rmse.mat) <- row.name
  
  if(remove == TRUE)
  {
    rmse.mat <- rmse.mat[,-4]
  }
  
  
  ##### sqrt_mse / min among different link functions #####
  rmse.comp <- t(apply(rmse.mat,1,function(x){x/min(x)}))
  rmse.comp <- cbind(index = 1:nrow(rmse.comp),rmse.comp)
  rmse.comp <- as.data.frame(rmse.comp)
  rmse.melt <- melt(rmse.comp,value.name = 'rmse',id.vars = 'index',variable.name = 'model')
  
  gp <- ggplot(data = rmse.melt,aes(x=index,y=rmse,color=model,linetype=model))+geom_line()+scale_x_discrete(labels=row.name) + theme_bw()+ylab('ratio') + xlab('assumption')
  
  res <- list(rmse = rmse.comp[,-1], gp=gp)
  return(res)
}

