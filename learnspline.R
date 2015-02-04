library(splines)
x <- rnorm(10)
b1 <- bs(x,knots= c(0.1,0.5),degree=1,intercept=TRUE)
b2 <- bs(x,knots=c(0.1,0.5),degree=2,intercept=TRUE)
kn <- c(rep(min(x),3),0.1,0.5,rep(max(x),3))

b1[,1]*(0.1-x)/(0.1-min(x)) == b2[,1]
b1[,1]*(x-min(x))/(0.1-min(x))+b1[,2]*(0.5-x)/(0.5-min(x)) == b2[,2]
b1[,2]*(x-min(x))/(0.5-min(x))+b1[,3]*(max(x)-x)/(max(x)-0.1)== b2[,3]
(x-0.1)*b1[,3]/(max(x)-0.1)+b1[,4]*(max(x)-x)/(max(x)-0.5) ==b2[,4]
(x-0.5)*b1[,4]/(max(x)-0.5) == b2[,5]


((0.1-x)>0)*(0.1-x)/(0.1-min(x))

((0.1-x)>0)*(x-min(x))/(0.1-min(x))+(x >0.1 & x<0.5)*(0.5-x)/(0.5-0.1)

b3 <- splineDesign(knots=c(min(x),min(x),min(x),0.1,0.5,max(x),max(x),max(x)), x, ord = 3, derivs=rep(0,10),outer.ok=TRUE)
b3.deriv <- splineDesign(knots=c(min(x),min(x),min(x),0.1,0.5,max(x),max(x),max(x)), x, ord = 3, derivs=rep(1,10),outer.ok=TRUE)


Mat <- matrix(c(-1,0,0,0,1,-1,0,0,0,1,-1,0,0,0,1,-1,0,0,0,1),nrow=4)
qvalue <- 2/c(0.1-min(x),0.5-min(x),max(x)-0.1,max(x)-0.5)
as.matrix(b1)%*%(qvalue*Mat%*%c(0,0,0,1,0))


#### consider degree =3 
bp2 <- bs(x,knots=c(0.1,0.5),degree=2,intercept=TRUE)
bp3 <- bs(x,knots=c(0.1,0.5),degree=3,intercept=TRUE)
b3.deriv <- splineDesign(knots=c(min(x),min(x),min(x),min(x),0.1,0.5,max(x),max(x),max(x),max(x)), x, ord = 4, derivs=rep(1,10),outer.ok=TRUE)


Mat <- matrix(c(-1,0,0,0,0,1,-1,0,0,0,0,1,-1,0,0,0,0,1,-1,0,0,0,0,1,-1,0,0,0,0,1),nrow=5)
qvalue <- 3/c(0.1-min(x),0.5-min(x),max(x)-min(x),max(x)-0.1,max(x)-0.5)
as.matrix(bp2)%*%(qvalue*Mat%*%c(1,0,0,0,0,0))
