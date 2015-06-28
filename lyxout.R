
library(xtable)
library(plyr)
create.bf <- function(x)
{
  x <- round(x,4)
  index <- which(rank(x) <= 3)
  
  x[index] <- paste0("\\textbf{", x[index], "}")
  return(x)
}

weights0 <- c('equal','both','left','right')

output.fun <- function(path,weight.arg,row.name)
{
  res.obj <- load(path)
  nw <- length(weight.arg)
  index <- which(weights0 %in% weight.arg)
 
  
  nres <- length(res.obj)
  
  prmse.ls <- ratio.ls <- list()
  
  for(w in 1:nw)
  {
    sub.fun <- function(ind)
    {
      colMeans(eval(parse(text=res.obj[ind]))$prmse.ls[[index[w]]])
    }
    
    prmse.mat <- ldply(1:nres,sub.fun)
    row.names(prmse.mat) <- row.name
    prmse.ls[[w]] <- prmse.mat

    
    ratio.ls[[w]] <- t(apply(prmse.mat,1,function(x){x/min(x)}))

   
  }
  
  names(prmse.ls) <-names(ratio.ls) <-  weight.arg
  
  return(list(prmse.ls = prmse.ls,ratio.ls = ratio.ls))
  
}


