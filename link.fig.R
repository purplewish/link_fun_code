###### consider splogit and gev #### distribution ####
setwd("/Users/Xin/Research/link_function/")
library(ggplot2)
library(gridExtra)
###gev ###
pdf('document/figures/comparison/gevlink.pdf',width = 8,height = 6)
label.gev <- c('logit',expression(paste(xi,'=',-1,sep=' ')),expression(paste(xi,'=',-0.3,sep=' ')),expression(paste(xi,'=',1,sep=' ')))
ggplot(data.frame(x=c(-7, 7)), aes(x)) + 
  stat_function(fun=function(x) plogis(x,location = 0,scale = 1),aes(color='logit'))+
  stat_function(fun = function(x) 1-pgev(-x,loc = 0,scale = 1,shape = -1), aes(color='gev1' ))+
  stat_function(fun = function(x) 1-pgev(-x,loc = 0,scale = 1,shape = -0.3), aes(color='gev2' ))+
  stat_function(fun = function(x) 1-pgev(-x,loc = 0,scale = 1,shape = 1), aes(color='gev3' ))+
  scale_colour_manual(name='link',values=c('logit'='black','gev1'="blue",'gev2'="green",'gev3'='red'),labels= label.gev,breaks=c('logit','gev1','gev2','gev3'))+theme_bw()+theme(legend.text.align=0)
dev.off()


splogit.link <- function(eta0,r0)
{
  if(r0 >0 & r0 <=1)
  {
    prob0 <- exp(eta0)/((1+exp(eta0/r0))^r0)
  }
  
  if(r0 > 1)
  {
    eta.new <- -r0*eta0
    prob0 <- 1-(exp(eta.new)/(1+exp(eta.new)))^(1/r0)
  }
  return(prob0)
}


#### splogit ####

pdf('document/figures/comparison/splogitlink.pdf',width = 8,height = 6)
label.splogit <- c(paste('r','=',0.1,sep=' '),paste('r','=',0.5,sep=' '),paste('r','=',2,sep=' '),paste('r','=','5',sep=''))

ggplot(data.frame(x=c(-7, 7)), aes(x)) + 
  stat_function(fun=function(x) plogis(x,location = 0,scale = 1),aes(color='logit'))+
  stat_function(fun = function(x) splogit.link(eta0 = x,r0 = 0.1), aes(color='splogit1'))+
  stat_function(fun = function(x) splogit.link(eta0 = x,r0 = 0.5), aes(color='splogit2'))+
  stat_function(fun = function(x) splogit.link(eta0 = x,r0 = 2), aes(color='splogit3'))+
  stat_function(fun = function(x) splogit.link(eta0 = x,r0 = 5), aes(color='splogit4'))+
  scale_colour_manual(name='link',values=c('logit'='black','splogit1'="blue",'splogit2'='green','splogit3'='red','splogit4'='orange'),labels=c('logit',label.splogit),breaks=c('logit',paste('splogit',1:4,sep='')))+theme_bw()+theme(legend.text.align=0)
dev.off()


### comparison between gev and splogit ####

r.fun <- function(r,xi,interval,nseq)
{
  x <- seq(interval[1],interval[2],length.out = nseq)
  gev.value <- 1-pgev(-x,loc=0,scale=1,shape=xi)
  splogit.value <- splogit.link(eta0 = x,r0 = r)
  dif <- sum((gev.value - splogit.value)^2)
  return(dif)
}

r1.fun <- function(r)
{
  r.fun(r,xi=-1,interval=c(-6,6),nseq=2000)
}

r2.fun <- function(r)
{
  r.fun(r,xi=1,interval=c(-5,5),nseq=2000)
}
optimize(r1.fun,c(0,10))
optimize(r2.fun,c(0.1,2))

r.fun.skew <- function(r,xi)
{
  if(r>0 & r <1)
  {
    value <- 1-1*(r/(r+1))^(r)-exp(-(1+xi))
  }
  
  if(r >1)
  {
    value<- r*(1/(r+1))^(1/r)- exp(-(1+xi))
  }
  return(value)
}

rs1 <- function(r)
{
  r.fun.skew(r,-1)
}


library(gridExtra)

pdf('document/figures/comparison/comparisonlink.pdf',width = 14,height = 6)
c1 <- ggplot(data.frame(x=c(-7, 7)), aes(x)) + 
  stat_function(fun = function(x) 1-pgev(-x,loc = 0,scale = 1,shape = -1), aes(color='gev2' ))+
  stat_function(fun = function(x) splogit.link(eta0 = x,r0 = 1.1174), aes(color='splogit2'))+
  stat_function(fun = function(x) splogit.link(eta0=x,r0=0.7840), aes(color='splogit1'))+
  scale_colour_manual(name='link',values=c('gev2'="green",'splogit2'='blue','splogit1'='red'),labels=c(expression(paste('GEV: ',xi,'=',-1)),paste('splogit:','r','=',1.1174),paste('splogit:','r','=',0.7840)))+
  theme_bw()+theme(legend.text.align=0)

c2 <- ggplot(data.frame(x=c(-7, 7)), aes(x)) + 
  stat_function(fun = function(x) 1-pgev(-x,loc = 0,scale = 1,shape = 1), aes(color='gev2' ))+
  stat_function(fun = function(x) splogit.link(eta0 = x,r0 = 0.4986), aes(color='splogit2'))+
  scale_colour_manual(name='link',values=c('gev2'="green",'splogit2'='blue'),labels=c(expression(paste('GEV: ',xi,'=',1)),paste('splogit:','r','=',0.4986)))+theme_bw()+
  theme(legend.text.align=0)
grid.arrange(c1,c2,ncol=2)
dev.off()



###### symmetric link #####
pdf("document/figures/symmetric_link.pdf",width = 6,height = 4)
label.sym <- c('logit',"probit",expression(paste(nu,'=',1,sep=' ')),expression(paste(nu,'=',2,sep=' ')))
ggplot(data.frame(x=c(-7, 7)), aes(x)) + 
  stat_function(fun=function(x) plogis(x,location = 0,scale = 1),aes(color='logit'))+
  stat_function(fun = function(x) pnorm(x,0,1),aes(color="probit"))+
  stat_function(fun = function(x) pt(x,df = 1),aes(color = "t1"))+
  stat_function(fun = function(x) pt(x,df = 2),aes(color = "t2"))+
  scale_colour_manual(name='link',values=c('logit'='black','probit'="blue",'t1'="green",'t2'='red'),labels= label.sym,breaks=c('logit','probit','t1','t2'))+theme_bw()+theme(legend.text.align=0)
dev.off()


pdf("document/figures/asymmetric_link.pdf",width = 6,height = 4)
label.asym <- c(paste('r','=',0.2,sep=' '),paste('r','=',5,sep=''),expression(paste(xi,'=',-0.5,sep=' ')),expression(paste(xi,'=',0.5,sep=' ')))
ggplot(data.frame(x=c(-7, 7)), aes(x)) + 
  stat_function(fun=function(x) plogis(x,location = 0,scale = 1),aes(color='logit'))+
  stat_function(fun = function(x) splogit.link(eta0 = x,r0 = 0.2), aes(color='splogit1'))+
  stat_function(fun = function(x) splogit.link(eta0 = x,r0 = 5), aes(color='splogit2'))+
  stat_function(fun = function(x) 1-pgev(-x,loc = 0,scale = 1,shape = -0.5), aes(color='gev1' ))+
  stat_function(fun = function(x) 1-pgev(-x,loc = 0,scale = 1,shape = 0.5), aes(color='gev2' ))+
  scale_colour_manual(name='link',values=c('logit'='black','splogit1'="blue",'splogit2'='green','gev1'='red','gev2'='orange'),labels=c('logit',label.asym),breaks=c('logit',paste('splogit',1:2,sep=''),paste("gev",1:2,sep="")))+theme_bw()+theme(legend.text.align=0)
dev.off()



##### curve ####

f1 <- function(x1)
{
  x1^2+x1-3
}

f2 <- function(x1)
{
  -0.2*(x1-3)^2-0.15*x1+3
}

f3 <- function(x1)
{
  0.2*(x1+1)^3-2
}


pdf("document/figures/curve.pdf",width = 6,height = 4)
ggplot(data.frame(x=c(-3, 3)), aes(x)) + 
  stat_function(fun = f1,aes(color='curve1'))+
  stat_function(fun = f2,aes(color='curve2'))+
  stat_function(fun = f3,aes(color='curve3'))+
  scale_colour_manual(name='curve',values=c('curve1'='black','curve2'="blue",'curve3'='green'),labels=c("curve1","curve2","curve3"),breaks=c("curve1","curve2","curve3"))+theme_bw()+theme(legend.text.align=0)
dev.off()

