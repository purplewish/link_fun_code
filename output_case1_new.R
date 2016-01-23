setwd("/Users/Xin/Research/link_function/")
source("/Users/Xin/Research/link_function/link_fun_code/lyxout.R") 
row.name <- c('logit','probit','robit(0.6)','robit(1)','robit(2)','gev(1)','gev(0.5)',
              'gev(-0.5)',"gev(-1)",'splogit(0.6)','splogit(1.5)')

weight.arg <- c("equal","both","left","right") 
truncated.arg <- c("greater","less","both")
pathn1_case1 <- 'output/deviance_gcv/output100_nonlinear_case1.RData' 
pathn2_case1 <- 'output/deviance_gcv/output200_nonlinear_case1.RData' 
pathn3_case1 <- 'output/deviance_gcv/output500_nonlinear_case1.RData'
resn100_case1 <- output.fun(path = pathn1_case1,weight.arg = weight.arg,
                            truncated.arg = truncated.arg,row.name = row.name) 
# resn200_case1 <- output.fun(path = pathn2_case2,weight.arg = weight.arg,row.name = row.name) 
# resn500_case1 <- output.fun(path = pathn3_case2,weight.arg = weight.arg,row.name = row.name)

###### 100####
print(xtable(t(apply(resn100_case1$ratio.ls$equal,1,create.bf)),caption = "sample size is 100 for RMSE of case 1",label = 'case1_100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn100_case1$ratio.ls$both,1,create.bf)),
             caption = "sample size is 100 for weighted RMSE of case 1 (1)",
             label = 'case1_100_1'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn100_case1$ratio.ls$data,1,create.bf)),
             caption = "sample size is 100 for weighted RMSE of case 1 (2)",
             label = 'case1_100_2'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn100_case1$ratio.ls$greater,1,create.bf)),
             caption = "sample size is 100 for weighted RMSE of case 1 (3)",
             label = 'case1_100_3'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)


######## plot ######
library(ggplot2)
pdf("document/figures/curve.pdf")
ggplot(data.frame(x=c(-3.5, 2.5)), aes(x))+
  stat_function(fun = function(x) x^2+x-3, aes(color='curve1'))+
  stat_function(fun=function(x)-0.2*(x-3)^2-0.15*x+3,aes(color='curve2'))+
  stat_function(fun=function(x)0.2*(x+1)^3-2,aes(color='curve3'))+
  scale_colour_manual(name='curve',values=c('curve1'="green",'curve2'='blue','curve3'='red'),labels=c('curve1','curve2','curve3'))+
  theme_bw()+theme(legend.text.align=0)
dev.off()