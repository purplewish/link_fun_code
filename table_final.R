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
pathn1 <- 'output/deviance_gcv/new/output100_binary_all.RData' 
pathn2 <- 'output/deviance_gcv/new/output200_binary_all.RData' 
pathn3 <- 'output/deviance_gcv/new/output500_binary_all.RData'
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



#### nonlinear case 1 and case 2

pathn1_case1 <- 'output/deviance_gcv/output100_nonlinear_case1_all.RData' 
pathn2_case1 <- 'output/deviance_gcv/output200_nonlinear_case1_all.RData' 
pathn3_case1 <- 'output/deviance_gcv/output500_nonlinear_case1_all.RData'

resn100_case1 <- output.fun(path = pathn1_case1,weight.arg = weight.arg,
                            truncated.arg = truncated.arg,row.name = row.namenon) 
resn200_case1 <- output.fun(path = pathn2_case1,weight.arg = weight.arg,
                            truncated.arg = truncated.arg,row.name = row.namenon) 
resn500_case1 <- output.fun(path = pathn3_case1,weight.arg = weight.arg,
                            truncated.arg = truncated.arg,row.name = row.namenon)

pathn1_case2 <- 'output/deviance_gcv/output100_nonlinear_case2_all.RData' 
pathn2_case2 <- 'output/deviance_gcv/output200_nonlinear_case2_all.RData' 
pathn3_case2 <- 'output/deviance_gcv/output500_nonlinear_case2_all.RData'

resn100_case2 <- output.fun(path = pathn1_case2,weight.arg = weight.arg,
                            truncated.arg = truncated.arg,row.name = row.namenon) 
resn200_case2 <- output.fun(path = pathn2_case2,weight.arg = weight.arg,
                            truncated.arg = truncated.arg,row.name = row.namenon) 
resn500_case2 <- output.fun(path = pathn3_case2,weight.arg = weight.arg,
                            truncated.arg = truncated.arg,row.name = row.namenon)

tab.final100non <- tab.final200non <- tab.final500non <-  data.frame(matrix(0,11,6))
weightuse <- c("equal","data","greater")
for(j in 1:3)
{
  tab.final100non[,j] <- first.three(resn100_case1$ratio.ls[[weightuse[j]]],shortnamenon)
  tab.final100non[,j+3] <- first.three(resn100_case2$ratio.ls[[weightuse[j]]],shortnamenon)
  
  tab.final200non[,j] <- first.three(resn200_case1$ratio.ls[[weightuse[j]]],shortnamenon)
  tab.final200non[,j+3] <- first.three(resn200_case2$ratio.ls[[weightuse[j]]],shortnamenon)
  
  tab.final500non[,j] <- first.three(resn500_case1$ratio.ls[[weightuse[j]]],shortnamenon)
  tab.final500non[,j+3] <- first.three(resn500_case2$ratio.ls[[weightuse[j]]],shortnamenon)
}

rownames(tab.final100non) <-  rownames(tab.final200non) <- rownames(tab.final500non) <-row.name
print(xtable(tab.final100non))
print(xtable(tab.final200non))
print(xtable(tab.final500non))



############ RMSE box plot ###########
library(reshape2)
library(ggplot2)
shortname1 <- c("L","P","R","G","S","PS","MPS","A")

row.namenon <- c('logit','probit','robit(0.6)','robit(1)','robit(2)','gev(1)','gev(0.5)','gev(-0.5)',"gev(-1)",'splogit(0.6)','splogit(1.5)')
shortnamenon1 <- c("L","P","R","G","S","PS","MPS")

res1 <- load(pathn1)

create.df.box <- function(path,row.name,shortnamenon1)
{
  obj <- load(path)
  nn <- length(row.name)
  df1 <- as.data.frame(matrix(0,100*nn,1+ length(shortnamenon1)))
  colnames(df1) <- c("Link",shortnamenon1)
  for(i in 1:nn)
  {
    index1 <- (i-1)*100 + 1
    index2 <- i*100
    df1[index1:index2,1] <- row.name[i]
    df1[index1:index2,-1] <- eval(parse(text=obj[i]))$prmse.ls$equal
  }
  
  dfm <- melt(df1,id.vars = "Link",variable.name = "assumption",value.name = "RMSE")
  dfm$Link <- factor(dfm$Link,levels = row.name)
  return(dfm)
}

df1 <- create.df.box(pathn1,row.name, shortname1)
df2 <- create.df.box(pathn2,row.name, shortname1)
df3 <- create.df.box(pathn3,row.name, shortname1)

pdf("document/STAT_COM/figures/linear100.pdf",width = 7.5,height = 10)
ggplot(df1,aes(x = assumption,y=RMSE))+
  geom_boxplot()+theme_bw() + theme(strip.background = element_blank())+
  facet_wrap(~Link,nrow = 4,ncol=3,scales = "free")
dev.off()

pdf("document/STAT_COM/figures/linear200.pdf",width = 7.5,height = 10)
ggplot(df2,aes(x = assumption,y=RMSE))+
  geom_boxplot()+theme_bw()+theme(strip.background = element_blank())+
  facet_wrap(~Link,nrow = 4,ncol=3, scales = "free")
dev.off()

pdf("document/STAT_COM/figures/linear500.pdf",width = 7.5,height = 10)
ggplot(df3,aes(x = assumption,y=RMSE))+
  geom_boxplot()+theme_bw()+theme(strip.background = element_blank())+
  facet_wrap(~Link,nrow = 4,ncol=3,scales = "free")
dev.off()


pdf("document/STAT_COM/figures/linear500_3fig.pdf",width = 7.5,height = 3)
ggplot(data = df3[df3$Link %in% c("logit","gev(-1)","splogit(0.2)"),],aes(x = assumption,y=RMSE))+
  geom_boxplot()+theme_bw()+theme(strip.background = element_blank())+
  facet_wrap(~Link,nrow = 1,ncol=3,scales = "free")
dev.off()

pathn3_case3 <- 'output/deviance_gcv/output500_nonlinear_case3_all.RData'
z1 <- load(pathn3_case3)

dfn3 <- create.df.box(pathn3_case3,row.namenon, shortnamenon1)

pdf("document/STAT_COM/figures/non_c3_500_3fig.pdf",width = 7.5,height = 3)
ggplot(data = dfn3[dfn3$Link %in% c("logit","gev(-1)","splogit(0.6)"),],aes(x = assumption,y=RMSE))+
  geom_boxplot()+theme_bw()+theme(strip.background = element_blank())+
  facet_wrap(~Link,nrow = 1,ncol=3,scales = "free")
dev.off()

