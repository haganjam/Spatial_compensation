
# miscellaneous functions

# taken from the rethinking package
# https://rdrr.io/github/rmcelreath/rethinking/src/R/utilities.r

# percentile confidence/credible interval
PCI <- function( samples , prob=0.89 ) {
  x <- sapply( prob , function(p) {
    a <- (1-p)/2
    quantile( samples , probs=c(a,1-a) )
  } )
  # now order inside-out in pairs
  n <- length(prob)
  result <- rep(0,n*2)
  for ( i in 1:n ) {
    low_idx <- n+1-i
    up_idx <- n+i
    # lower
    result[low_idx] <- x[1,i]
    # upper
    result[up_idx] <- x[2,i]
    # add names
    a <- (1-prob[i])/2
    names(result)[low_idx] <- paste(round(a*100,0),"%")
    names(result)[up_idx] <- paste(round((1-a)*100,0),"%")
  }
  return(result)
}
PI <- PCI

### END
