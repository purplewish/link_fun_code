setwd("/Users/Xin/Research/link_function/")
source("/Users/Xin/Research/link_function/link_fun_code/lyxout.R") 
row.name <- c('logit','probit','robit(0.6)','robit(1)','robit(2)','gev(1)','gev(0.5)',
              'gev(-0.5)',"gev(-1)",'splogit(0.6)','splogit(1.5)')

weight.arg <- c("equal","both","left","right","data") 
truncated.arg <- c("greater","less","both")
pathn1_case1 <- 'output/deviance_gcv/output100_nonlinear_case1_all.RData' 
pathn2_case1 <- 'output/deviance_gcv/output200_nonlinear_case1_all.RData' 
pathn3_case1 <- 'output/deviance_gcv/output500_nonlinear_case1_all.RData'
resn100_case1 <- output.fun(path = pathn1_case1,weight.arg = weight.arg,
                            truncated.arg = truncated.arg,row.name = row.name) 

resn200_case1 <- output.fun(path = pathn2_case1,weight.arg = weight.arg,
                            truncated.arg = truncated.arg,row.name = row.name)  
resn500_case1 <- output.fun(path = pathn3_case1,weight.arg = weight.arg,
                            truncated.arg = truncated.arg,row.name = row.name) 
###### 100####
print(xtable(t(apply(resn100_case1$ratio.ls$equal,1,create.bf)),caption = "RMSE of nonlinear case 1 with sample size 100",label = 'case1_100'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn100_case1$ratio.ls$data,1,create.bf)),
             caption = "Weighted RMSE of nonlinear case 1 with sample size 100",
             label = 'case1_100_2'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn100_case1$ratio.ls$greater,1,create.bf)),
             caption = "Absolute error of nonlinear case 1 with sample size 100",
             label = 'case1_100_3'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

## 200
print(xtable(t(apply(resn200_case1$ratio.ls$equal,1,create.bf)),caption = "RMSE of nonlinear case 1 with sample size 200",label = 'case1_200_1'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn200_case1$ratio.ls$data,1,create.bf)),
             caption = "Weighted RMSE of nonlinear case 1 with sample size 200",
             label = 'case1_200_2'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn200_case1$ratio.ls$greater,1,create.bf)),
             caption = "Absolute error of nonlinear case 1 with sample size 200",
             label = 'case1_200_3'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

## 500
print(xtable(t(apply(resn500_case1$ratio.ls$equal,1,create.bf)),caption = "RMSE of nonlinear case 1 with sample size 500",label = 'case1_500_1'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn500_case1$ratio.ls$data,1,create.bf)),
             caption = "Weighted RMSE of nonlinear case 1 with sample size 500",
             label = 'case1_500_2'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

print(xtable(t(apply(resn500_case1$ratio.ls$greater,1,create.bf)),
             caption = "Absolute error of nonlinear case 1 with sample size 500",
             label = 'case1_500_3'),caption.placement='top',table.placement = 'H', sanitize.text.function = function(x) x)

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