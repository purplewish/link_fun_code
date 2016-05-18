setwd("/Users/Xin/Research/link_function/")
source("/Users/Xin/Research/link_function/link_fun_code/lyxout.R") 
row.name <- c('logit','probit','robit(0.6)','robit(1)','robit(2)','gev(1)','gev(0.5)',
              'gev(-0.5)',"gev(-1)",'splogit(0.2)','splogit(5)')

weight.arg <- c("equal","both","left","right","data")
truncated.arg <- c("greater","less","both")
pathn1 <- 'output/deviance_gcv/new/output100_binary_all.RData' 
pathn2 <- 'output/deviance_gcv/new/output200_binary_all.RData' 
pathn3 <- 'output/deviance_gcv/new/output500_binary_all.RData'
resn100 <- output.fun(path = pathn1,weight.arg = weight.arg,
                      truncated.arg = truncated.arg,row.name = row.name) 
resn200 <- output.fun(path = pathn2,weight.arg = weight.arg,
                      truncated.arg = truncated.arg,row.name = row.name) 
resn500<- output.fun(path = pathn3,weight.arg = weight.arg,
                     truncated.arg = truncated.arg,row.name = row.name)

###### 100####
print(xtable(t(apply(resn100$ratio.ls$equal,1,create.bf)),caption = "sample size is 100 for RMSE of linear case",label = 'binary_100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn100$ratio.ls$both,1,create.bf)),
             caption = "sample size is 100 for weighted RMSE of case 1 (1)",
             label = 'binary_100_1'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn100$ratio.ls$data,1,create.bf)),
             caption = "sample size is 100 for weighted RMSE of case 1 (2)",
             label = 'binary_100_2'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn100$ratio.ls$greater,1,create.bf)),
             caption = "sample size is 100 for weighted RMSE of case 1 (3)",
             label = 'binary_100_3'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

### 200
print(xtable(t(apply(resn200$ratio.ls$equal,1,create.bf)),caption = "sample size is 200 for RMSE of linear case",label = 'binary_200'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn200$ratio.ls$both,1,create.bf)),
             caption = "sample size is 200 for weighted RMSE of case 1 (1)",
             label = 'binary_200_1'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn200$ratio.ls$data,1,create.bf)),
             caption = "sample size is 200 for weighted RMSE of case 1 (2)",
             label = 'binary_200_2'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn200$ratio.ls$greater,1,create.bf)),
             caption = "sample size is 200 for weighted RMSE of case 1 (3)",
             label = 'binary_200_3'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)



###### 500####
print(xtable(t(apply(resn500$ratio.ls$equal,1,create.bf)),caption = "sample size is 500 for RMSE of linear case",label = 'binary_500'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn500$ratio.ls$both,1,create.bf)),
             caption = "sample size is 500 for weighted RMSE of case 1 (1)",
             label = 'binary_500_1'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn500$ratio.ls$data,1,create.bf)),
             caption = "sample size is 500 for weighted RMSE of case 1 (2)",
             label = 'binary_500_2'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn500$ratio.ls$greater,1,create.bf)),
             caption = "sample size is 500 for weighted RMSE of case 1 (3)",
             label = 'binary_500_3'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)


