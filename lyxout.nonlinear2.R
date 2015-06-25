
create.bf <- function(x)
{
  x <- round(x,4)
  index <- which(rank(x) <= 3)
  
  x[index] <- paste0("\\textbf{", x[index], "}")
  return(x)
}
source('link_fun_code/tab.fig.fun.R')
col.name <-  c('logit','probit','robit','gev','splogit','pspline') 
row.name <- c('logit','probit','robit(0.6)','robit(1)','robit(2)','gev(1)','gev(0.5)','gev(-0.5)',"gev(-1)",'splogit(0.6)','splogit(1.5)')

load('output/output100_nonlinear2.RData')
resp <- tab.fig.fun(prmse.out,col.name=col.name,row.name=row.name,remove=FALSE)
res <- apply(resp$rmse,1,create.bf)
print(xtable(t(res),caption = "sample size is 100 for RMSE in nonlinear case 2",label = 'nns100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)



load('output/output200_nonlinear2.RData')
resp <- tab.fig.fun(prmse.out,col.name=col.name,row.name=row.name,remove=FALSE)
res <- apply(resp$rmse,1,create.bf)
print(xtable(t(res),caption = "sample size is 200 for RMSE in nonlinear case 2",label = 'nns200'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)


load('output/output500_nonlinear2.RData')
resp <- tab.fig.fun(prmse.out,col.name=col.name,row.name=row.name,remove=FALSE)
res <- apply(resp$rmse,1,create.bf)
print(xtable(t(res),caption = "sample size is 500 for RMSE in nonlinear case 2",label = 'nns500'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)



