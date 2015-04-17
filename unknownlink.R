#### algorithm #### 
#### something wrong with identification ####
library(actuar)
deg <- 3
nknots <- 6
tol <- 1e-8
monotone = TRUE
boundary = c(-5,5)

### standardized x and determine the value of a###
xmat <- cbind(x1,x2)
xmats <- scale(xmat)

atu <- quantile(sqrt(rowSums(xmats^2)),0.999)

### initial value of beta ####
glm.fit.logit <- glm(y0~0+x1+x2,family = binomial())
beta00 <- coef(glm.fit.logit)
beta00 <- beta00/sqrt(sum(beta00^2))
eta.old <- xmats%*%beta00
q.old <- (eta.old/atu+1)/2
d.value <- 2
Ut <- pgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )

nknots <- 7
bs.nc <- nknots+deg-1
delta0 <- rep(1/bs.nc,bs.nc)
knots <- seq(0,1,length.out = nknots)

beta.old <- beta00
delta.old <- delta0
link.est <- function(y0,x0,deg = 3,nknots,monotone=TRUE,beta0,delta0,tol = 1e-8,boundary)
{
  xmat <- cbind(1,x0)
  beta.old <- beta0
  delta.old <- delta0

 ##### link function #### 
  repeat
    {   
      beta.comp <- beta.old
      delta.comp <- delta.old
      
      eta.old <- xmats%*%beta.comp
      q.old <- (eta.old/atu+1)/2
      d.value <- 2
      Ut <- pgenbeta(q.old,shape1 = (d.value+1)/2,shape2 = (d.value+1)/2,shape3 = 1,scale = 1  )
      bs0 <- bs(Ut,knots=knots[c(-1,-length(knots))],degree=deg,intercept=TRUE)
      delta.old <- coef(glm(y0~0+bs0,family = 'binomial'))
      
    if(!monotone)
    {
      repeat
      {       
        bs.eta <- bs0%*%delta.old
        bs.mu <- exp(bs.eta)/(1+exp(bs.eta))
        wt <- diag(as.numeric(bs.mu*(1-bs.mu)))
        z <- bs.eta+(y0-bs.mu)/as.numeric(bs.mu*(1-bs.mu)) 
        delta.update <- solve(t(bs0)%*%wt%*%bs0)%*%t(bs0)%*%wt%*%z
        print(delta.update)
        delta.old <- delta.update
        if(diff.value <= tol){break}
      }  
    }
    
    if(monotone)
    {
      kp <- 1e6
      Dmat <- matrix(0,ncol=bs.nc,nrow=bs.nc-1)
      diag(Dmat) <- -1
      Dmat[cbind(1:(bs.nc-1),2:bs.nc)] <- 1
      
      repeat
      {
        bs.eta <- bs0%*%delta.old
        bs.mu <- exp(bs.eta)/(1+exp(bs.eta))
        wt <- diag(as.numeric(bs.mu*(1-bs.mu)))
        z <- bs.eta+(y0-bs.mu)/as.numeric(bs.mu*(1-bs.mu)) 
        Vmat <- diag(as.numeric(Dmat%*%delta.old) <0 ) 
        delta.update <- ginv(t(bs0)%*%wt%*%bs0 + kp*t(Dmat)%*%Vmat%*%Dmat)%*%t(bs0)%*%wt%*%z
        diff.value <- mean((delta.update-delta.old)^2)
        print(delta.update)
        delta.old <- delta.update
        if(diff.value <= tol){break}
        
      }
      
    }
    
    ### update beta ####
    
    repeat
    {
      eta.old<- xmat%*%beta.old
      #eta.stand <- (eta.old - mean(eta.old))/sqrt(mean((eta.old - mean(eta.old))^2))
      eta.stand <- eta.old
      bs.old <- bs(eta.stand,knots=knots[c(-1,-length(knots))],degree=deg,Boundary.knots =c(-4,4),intercept=TRUE)
      muhat <- exp(bs.old %*% delta.update)/(1+exp(bs.old %*% delta.update))
      bs.deriv <- splineDesign(knots= c(rep(-4,4),knots[c(-1,-length(knots))],rep(4,4)), eta.stand, ord = 4, derivs=rep(1,length(x0)),outer.ok=TRUE)
      fun.deriv <- bs.deriv%*% delta.update
      du.eta <- muhat*(1-muhat)* fun.deriv
      z.beta <- eta.stand +(y0 - muhat)/du.eta
      wt.beta <- diag(as.numeric(muhat*(1-muhat)*fun.deriv^2))
      beta.update <- solve(t(xmat)%*%wt.beta%*%xmat)%*%t(xmat)%*%wt.beta%*%z.beta
      diff.value <- mean((beta.update-beta.old)^2)
      beta.old <- beta.update
      if(diff.value <= tol){break} 
    }
    
    diff.beta <- mean()    
  }

}




#### final result 
plot(x0,prob0,type='l')
lines(x0,muhat,col='red')  
