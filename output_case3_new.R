setwd("C:/Users/xinwang/Research/link_function/")
source("C:/Users/xinwang/Research/link_function/link_fun_code/lyxout.R") 
row.name <- c('logit','probit','robit(0.6)','robit(1)','robit(2)','gev(1)','gev(0.5)',
              'gev(-0.5)',"gev(-1)",'splogit(0.6)','splogit(1.5)')

weight.arg <- c("equal","both","left","right") 
pathn1_case3 <- 'output/new/output100_nonlinear_case3.RData' 
pathn2_case3 <- 'output/new/output200_nonlinear_case3.RData' 
pathn3_case3 <- 'output/new/output500_nonlinear_case3.RData'
resn100_case3 <- output.fun(path = pathn1_case3,weight.arg = weight.arg,row.name = row.name) 
resn200_case3 <- output.fun(path = pathn2_case3,weight.arg = weight.arg,row.name = row.name) 
resn500_case3 <- output.fun(path = pathn3_case3,weight.arg = weight.arg,row.name = row.name)

###### 100####
print(xtable(t(apply(resn100_case3$ratio.ls$equal,1,create.bf)),caption = "sample size is 100 for RMSE of case 3",
             label = 'case3_100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn100_case3$ratio.ls$both,1,create.bf)),
             caption = "sample size is 100 for weighted RMSE of case 3 (1)",
             label = 'case3_100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn100_case3$ratio.ls$left,1,create.bf)),
             caption = "sample size is 100 for weighted RMSE of case 3 (2)",
             label = 'case3_100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn100_case3$ratio.ls$right,1,create.bf)),
             caption = "sample size is 100 for weighted RMSE of case 3 (3)",
             label = 'case3_100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

###### 200####
print(xtable(t(apply(resn200_case3$ratio.ls$equal,1,create.bf)),caption = "sample size is 200 for RMSE of case 3",
             label = 'case3_200'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn200_case3$ratio.ls$both,1,create.bf)),
             caption = "sample size is 200 for weighted RMSE of case 3 (1)",
             label = 'case3_200'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn200_case3$ratio.ls$left,1,create.bf)),
             caption = "sample size is 200 for weighted RMSE of case 3 (2)",
             label = 'case3_200'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn200_case3$ratio.ls$right,1,create.bf)),
             caption = "sample size is 200 for weighted RMSE of case 3 (3)",
             label = 'case3_200'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

###### 500####
print(xtable(t(apply(resn500_case3$ratio.ls$equal,1,create.bf)),caption = "sample size is 500 for RMSE of case 3",
             label = 'case3_500'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn500_case3$ratio.ls$both,1,create.bf)),
             caption = "sample size is 500 for weighted RMSE of case 3 (1)",
             label = 'case3_500'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn500_case3$ratio.ls$left,1,create.bf)),
             caption = "sample size is 500 for weighted RMSE of case 3 (2)",
             label = 'case3_500'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn500_case3$ratio.ls$right,1,create.bf)),
             caption = "sample size is 500 for weighted RMSE of case 3 (3)",
             label = 'case3_500'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)