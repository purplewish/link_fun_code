setwd("C:/Users/xinwang/Research/link_function/")
source("C:/Users/xinwang/Research/link_function/link_fun_code/lyxout.R") 
row.name <- c('logit','probit','robit(1)','robit(2)','robit(0.6)','gev(0.5)',
              'gev(-0.5)',"gev(-1)",'splogit(0.6)','splogit(1.5)')

weight.arg <- c("equal","both","left","right") 
path_paper<- 'output/output100_binary_paper.RData' 

res_paper <- output.fun(path = path_paper,weight.arg = weight.arg,row.name = row.name) 
# resn200_case1 <- output.fun(path = pathn2_case2,weight.arg = weight.arg,row.name = row.name) 
# resn500_case1 <- output.fun(path = pathn3_case2,weight.arg = weight.arg,row.name = row.name)

###### 100####
print(xtable(t(apply(res_paper$ratio.ls$equal,1,create.bf)),caption = "sample size is 100 for RMSE",label = 'paper_100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(res_paper$ratio.ls$both,1,create.bf)),caption = "sample size is 100 for weighted RMSE (1) ",label = 'paper_100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(res_paper$ratio.ls$left,1,create.bf)),caption = "sample size is 100 for weighted RMSE (2)",label = 'paper_100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(res_paper$ratio.ls$right,1,create.bf)),caption = "sample size is 100 for weighted RMSE (3)",label = 'paper_100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

