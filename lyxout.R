
library(xtable)
library(plyr)
create.bf <- function(x)
{
  x <- round(x,2)
  index <- which(rank(x) <= 3)
  
  x[index] <- paste0("\\textbf{", x[index], "}")
  return(x)
}

weights0 <- c('equal','both','left','right',"data")
truncated0 <- c('greater',"less","both")

output.fun <- function(path,weight.arg,truncated.arg,row.name)
{
  res.obj <- load(path)
  nw <- length(weight.arg)
  nt <- length(truncated.arg)
  index <- which(weights0 %in% weight.arg)
  indext <- which(truncated0 %in% truncated.arg)
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
  
  for(w in 1:nt)
  {
    sub.funt <- function(ind)
    {
      colMeans(eval(parse(text=res.obj[ind]))$truncated.ls[[indext[w]]])
    }
    prmse.mat <- ldply(1:nres,sub.funt)
    row.names(prmse.mat) <- row.name
    prmse.ls[[w+nw]] <- prmse.mat
    ratio.ls[[w+nw]] <- t(apply(prmse.mat,1,function(x){x/min(x)}))
  }
  
  names(prmse.ls) <-names(ratio.ls) <-  c(weight.arg,truncated.arg)
  
  return(list(prmse.ls = prmse.ls,ratio.ls = ratio.ls))
  
}


first.three <- function(res.ratio,shortname)
{
  res.ratio <- round(res.ratio,4)
  sub.fun <- function(x)
  {
    rankx <- rank(x)
    threemodel <- shortname[sort(rankx,index.return=TRUE)$ix][1:3]
    res <- paste(threemodel,collapse="/")
    return(res)
  }
  res0 <- apply(res.ratio,1,sub.fun)
  return(res0)
}
