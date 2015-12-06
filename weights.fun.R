weights.fun <- function(arg,data)
{
  minx <- min(data)
  maxx <- max(data)
  if(arg == 'equal')
  {
    weights <- 1/length(data)
  }
  if(arg == 'both')
  {
    weights = 1/dnorm(data,mean = -0.5)
    weights <- weights/sum(weights)
  }
  
  if(arg == 'left')
  {
    
    weights <- dchisq(data-minx,2)
    weights <- weights/sum(weights)
  }
  
  if(arg == 'right')
  {
    weights <- dchisq(data-minx,9)
    weights <- weights/sum(weights)
  }
  
  if(arg == "data")
  {
    weights <- dnorm(data,mean=-0.5)
    weights <- weights/sum(weights)
  }
  
  return(weights)
}