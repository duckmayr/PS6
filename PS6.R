sg.int <- function(g, lower, upper, dimensions, parallel=TRUE, ...){
  require("SparseGrid")
  
  lower <- floor(lower)
  upper <- ceiling(upper)
  
  if ( any(lower > upper) ) {
    stop("lower must be smaller than upper")
  }
  
  gridss <- as.matrix(expand.grid(lapply(1:dimensions, function(x){
    seq(lower[x], upper[x]-1, by=1)
  })))
  sp.grid <- createIntegrationGrid( 'KPU', dimension=dimensions, k=5 )
  nodes <- gridss[1,] + sp.grid$nodes
  weights <- sp.grid$weights
  
  for ( i in 2:nrow(gridss) ) {
    nodes <- rbind(nodes, gridss[i,] + sp.grid$nodes)  
    weights <- c(weights, sp.grid$weights)
  }
  
  if ( parallel ) {
    require('parallel')
    clus <- makeCluster(max(detectCores() - 1, 1))
    
    gx.sp <- parApply(cl=clus, X=nodes, MARGIN=1, FUN=g, ...)
    val.sp <- gx.sp %*% weights
    
    stopCluster(clus)
  } else {
    gx.sp <- apply(X=nodes, MARGIN=1, FUN=g, ...)
    val.sp <- gx.sp %*% weights
  }
  
  return(val.sp)
}

## Tests:

library(mvtnorm)
library(testthat)

myNorm <- function(x, dimensions){
  dmvnorm(x, mean=rep(0, dimensions), sigma=diag(rep(1, dimensions)))
}

realAns <- function(lower, upper, dimensions){
  return(as.numeric(pmvnorm(rep(lower, dimensions), rep(upper, dimensions),
                            rep(0, dimensions), diag(rep(1, dimensions)))))
}

comparison <- function(lower, upper, dimensions){
  return(as.numeric(sg.int(myNorm, rep(lower, dimensions),
                           rep(upper, dimensions), dimensions,
                           FALSE, dimensions)))
}

expect_equivalent(realAns(-2, 2, 2), comparison(-2, 2, 2))
expect_equivalent(realAns(-2, 2, 4), comparison(-2, 2, 4))

mySine <- function(x){
  return(sum(sin(x)))
}

mcInt <- function(a, b, n, dimensions){
  u <- runif(n*dimensions, a, b)
  these <- matrix(u, ncol=dimensions)
  x <- sapply(1:nrow(these), function(x) mySine(c(these[x, ])))
  return(mean(x)*(b-a)^dimensions)
}

comparison <- function(lower, upper, dimensions){
  return(as.numeric(sg.int(mySine, lower=rep(lower, dimensions),
                           upper=rep(upper, dimensions),
                           dimensions=dimensions, parallel=FALSE)))
}

expect_equal(mcInt(0, 9, 100000, 2), comparison(0, 9, 2), tol=0.1)
expect_equal(mcInt(0, 9, 100000, 3), comparison(0, 9, 3), tol=0.1)

expect_error(sg.int(mySine, lower='rep(0, 3)', upper=rep(9, 3), dimensions=3, parallel=FALSE),
             'non-numeric argument to mathematical function')
expect_error(sg.int(mySine, lower=rep(0, 3), upper=rep(9, 3), dimensions='3', parallel=FALSE),
             "wrong sign in 'by' argument")
expect_error(sg.int(mySine, lower=rep(9, 3), upper=rep(0, 3), dimensions=3, parallel=FALSE),
             "lower must be smaller than upper")
