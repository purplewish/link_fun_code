
library(xtable)
create.bf <- function(x)
{
  x <- round(x,4)
  index <- which(rank(x) <= 3)
  
  x[index] <- paste0("\\textbf{", x[index], "}")
  return(x)
}
source('link_fun_code/tab.fig.fun.R')
col.name <-  c('logit','probit','robit','gev','splogit','gam','pspline')
row.name <- c('logit','probit','robit(0.6)','robit(1)','robit(2)','gev(1)','gev(0.5)','gev(-0.5)',"gev(-1)",'splogit(0.2)','splogit(5)')

load('output/output100_binarynew.RData')



output.fun <- function(path,weight.arg,linear=TRUE,col.name,row.name)
{
  res.obj <- load(path)
  nw <- length(weight.arg)
  weights.use <- weights0[weights0 %in% weight.arg]
  
  nres <- length(res.obj)
  
  eval(parse(text=res.obj[1]))$prmse.ls[[2]]
  
}

for(w in 1:4)
{
  prmse.mat <- out.logit$prmse.ls[[w]],
}

resp <- tab.fig.fun(prmse.out,col.name=col.name,row.name=row.name,remove=FALSE)
res <- apply(resp$rmse,1,create.bf)
print(xtable(t(res),caption = "sample size is 100 for RMSE",label = 'ns100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

resw <- tab.fig.fun(wprmse.out,col.name=col.name,row.name=row.name,remove=FALSE)
res <- apply(resw$rmse,1,create.bf)
print(xtable(t(res),caption = "sample size is 100 for weighted RMSE",label = 'ns100w'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)


load('output/output200_binary.RData')
resp <- tab.fig.fun(prmse.out,col.name=col.name,row.name=row.name,remove=FALSE)
res <- apply(resp$rmse,1,create.bf)
print(xtable(t(res),caption = "sample size is 200 for RMSE",label = 'ns200'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

resw <- tab.fig.fun(wprmse.out,col.name=col.name,row.name=row.name,remove=FALSE)
res <- apply(resw$rmse,1,create.bf)
print(xtable(t(res),caption = "sample size is 200 for weighted RMSE",label = 'ns200w'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

load('output/output500_binary.RData')
resp <- tab.fig.fun(prmse.out,col.name=col.name,row.name=row.name,remove=FALSE)
res <- apply(resp$rmse,1,create.bf)
print(xtable(t(res),caption = "sample size is 500 for RMSE",label = 'ns500'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

resw <- tab.fig.fun(wprmse.out,col.name=col.name,row.name=row.name,remove=FALSE)
res <- apply(resw$rmse,1,create.bf)
print(xtable(t(res),caption = "sample size is 500 for weighted RMSE",label = 'ns500w'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)



