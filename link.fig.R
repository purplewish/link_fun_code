###### consider splogit and gev #### distribution ####
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
pdf('document/figures/comparison/comparisonlink.pdf',width = 14,height = 6)
c1 <- ggplot(data.frame(x=c(-7, 7)), aes(x)) + 
  stat_function(fun = function(x) 1-pgev(-x,loc = 0,scale = 1,shape = -1), aes(color='gev2' ))+
  stat_function(fun = function(x) splogit.link(eta0 = x,r0 = 2), aes(color='splogit2'))+
  scale_colour_manual(name='link',values=c('gev2'="green",'splogit2'='blue'),labels=c(expression(paste('GEV: ',xi,'=',-1)),paste('splogit:','r','=',2)))+
  theme_bw()+theme(legend.text.align=0)

c2 <- ggplot(data.frame(x=c(-7, 7)), aes(x)) + 
  stat_function(fun = function(x) 1-pgev(-x,loc = 0,scale = 1,shape = 1), aes(color='gev2' ))+
  stat_function(fun = function(x) splogit.link(eta0 = x,r0 = 0.02), aes(color='splogit2'))+
  scale_colour_manual(name='link',values=c('gev2'="green",'splogit2'='blue'),labels=c(expression(paste('GEV: ',xi,'=',1)),paste('splogit:','r','=',0.02)))+theme_bw()+
  theme(legend.text.align=0)
grid.arrange(c1,c2,ncol=2)
dev.off()