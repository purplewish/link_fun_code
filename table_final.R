setwd("/Users/Xin/Research/link_function/")
source("/Users/Xin/Research/link_function/link_fun_code/lyxout.R") 
row.name <- c('logit','probit','robit(0.6)','robit(1)','robit(2)','gev(1)','gev(0.5)',
              'gev(-0.5)',"gev(-1)",'splogit(0.2)','splogit(5)')
shortname <- c("L","P","R","G","S","PW","PM","A")
row.namenon <- c('logit','probit','robit(0.6)','robit(1)','robit(2)','gev(1)','gev(0.5)',
              'gev(-0.5)',"gev(-1)",'splogit(0.6)','splogit(1.5)')
shortnamenon <- c("L","P","R","G","S","PW","PM")

weight.arg <- c("equal","both","left","right","data")
truncated.arg <- c("greater","less","both")
pathn1 <- 'output/deviance_gcv/output100_binary_all.RData' 
pathn2 <- 'output/deviance_gcv/output200_binary_all.RData' 
pathn3 <- 'output/deviance_gcv/output500_binary_all.RData'
resn100 <- output.fun(path = pathn1,weight.arg = weight.arg,
                      truncated.arg = truncated.arg,row.name = row.name) 
resn200 <- output.fun(path = pathn2,weight.arg = weight.arg,
                      truncated.arg = truncated.arg,row.name = row.name) 
resn500<- output.fun(path = pathn3,weight.arg = weight.arg,
                     truncated.arg = truncated.arg,row.name = row.name)

pathn1_case3 <- 'output/deviance_gcv/output100_nonlinear_case3_all.RData' 
pathn2_case3 <- 'output/deviance_gcv/output200_nonlinear_case3_all.RData' 
pathn3_case3 <- 'output/deviance_gcv/output500_nonlinear_case3_all.RData'

resn100_case3 <- output.fun(path = pathn1_case3,weight.arg = weight.arg,
                            truncated.arg = truncated.arg,row.name = row.namenon) 
resn200_case3 <- output.fun(path = pathn2_case3,weight.arg = weight.arg,
                            truncated.arg = truncated.arg,row.name = row.namenon) 
resn500_case3 <- output.fun(path = pathn3_case3,weight.arg = weight.arg,
                            truncated.arg = truncated.arg,row.name = row.namenon)
tab.final100 <- tab.final200 <- tab.final500 <-  data.frame(matrix(0,11,6))
weightuse <- c("equal","data","greater")
for(j in 1:3)
{
  tab.final100[,j] <- first.three(resn100$ratio.ls[[weightuse[j]]],shortname)
  tab.final100[,j+3] <- first.three(resn100_case3$ratio.ls[[weightuse[j]]],shortnamenon)
  
  tab.final200[,j] <- first.three(resn200$ratio.ls[[weightuse[j]]],shortname)
  tab.final200[,j+3] <- first.three(resn200_case3$ratio.ls[[weightuse[j]]],shortnamenon)
  
  tab.final500[,j] <- first.three(resn500$ratio.ls[[weightuse[j]]],shortname)
  tab.final500[,j+3] <- first.three(resn500_case3$ratio.ls[[weightuse[j]]],shortnamenon)
}

rownames(tab.final100) <-  rownames(tab.final200) <- rownames(tab.final500) <-row.name
print(xtable(tab.final100))
print(xtable(tab.final200))
print(xtable(tab.final500))

