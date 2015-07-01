setwd("C:/Users/xinwang/Research/link_function/")
source("C:/Users/xinwang/Research/link_function/link_fun_code/lyxout.R") 
row.name <- c('logit','probit','robit(1)','robit(2)','robit(0.6)',
              'gev(1)','gev(0.5)','gev(-0.5)',"gev(-1)",'splogit(0.2)','splogit(5)') 
weight.arg <- c("equal","both","left","right") 
path1 <- 'output/output100_binarynew.RData' 
path2 <- 'output/output200_binarynew.RData' 
path3 <- 'output/output500_binarynew.RData'
res100 <- output.fun(path = path1,weight.arg = weight.arg,row.name = row.name) 
res200 <- output.fun(path = path2,weight.arg = weight.arg,row.name = row.name) 
res500 <- output.fun(path = path3,weight.arg = weight.arg,row.name = row.name)
print(xtable(t(apply(res100$ratio.ls$equal,1,create.bf)),caption = "sample size is 100 for RMSE",label = 'ns100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(res100$ratio.ls$both,1,create.bf)),caption = "sample size is 100 for weighted RMSE (1)",label = 'ns100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(res100$ratio.ls$left,1,create.bf)),caption = "sample size is 100 for weighted RMSE (2)",label = 'ns100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(res100$ratio.ls$right,1,create.bf)),caption = "sample size is 100 for weighted RMSE (3)",label = 'ns100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

##### 200####
print(xtable(t(apply(res200$ratio.ls$equal,1,create.bf)),caption = "sample size is 200 for RMSE",label = 'ns100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(res200$ratio.ls$both,1,create.bf)),caption = "sample size is 200 for weighted RMSE (1)",label = 'ns100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(res200$ratio.ls$left,1,create.bf)),caption = "sample size is 200 for weighted RMSE (2)",label = 'ns100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(res200$ratio.ls$right,1,create.bf)),caption = "sample size is 200 for weighted RMSE (3)",label = 'ns100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)


##### 500 #####
print(xtable(t(apply(res500$ratio.ls$equal,1,create.bf)),caption = "sample size is 500 for RMSE",label = 'ns100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(res500$ratio.ls$both,1,create.bf)),caption = "sample size is 500 for weighted RMSE (1)",label = 'ns100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(res500$ratio.ls$left,1,create.bf)),caption = "sample size is 500 for weighted RMSE (2)",label = 'ns100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(res500$ratio.ls$right,1,create.bf)),caption = "sample size is 500 for weighted RMSE (3)", label = 'ns100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)