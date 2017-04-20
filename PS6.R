sg.int <- function(g, lower, upper, dimensions, parallel=TRUE, ...){
  require("SparseGrid")
  
  # Make sure the upper and lower limits of integration are proper:
  lower <- floor(lower)
  upper <- ceiling(upper)
  if ( any(lower > upper) ) {
    stop("lower must be smaller than upper")
  }
  
  # This creates a matrix of all possible combinations of variable values,
  # like a Cartesian product
  gridss <- as.matrix(expand.grid(lapply(1:dimensions, function(x){
    seq(lower[x], upper[x]-1, by=1)
  })))
  
  # The next seven lines create nodes and weights for integration
  sp.grid <- createIntegrationGrid( 'KPU', dimension=dimensions, k=5 )
  nodes <- gridss[1,] + sp.grid$nodes
  weights <- sp.grid$weights
  for ( i in 2:nrow(gridss) ) {
    nodes <- rbind(nodes, gridss[i,] + sp.grid$nodes)  
    weights <- c(weights, sp.grid$weights)
  }
  
  # If the user wants to do parallel calculation, we do that
  if ( parallel ) {
    require('parallel')
    clus <- makeCluster(max(detectCores() - 1, 1))
    
    gx.sp <- parApply(cl=clus, X=nodes, MARGIN=1, FUN=g, ...)
    val.sp <- gx.sp %*% weights
    
    stopCluster(clus)
  } else { # If not, we don't
    gx.sp <- apply(X=nodes, MARGIN=1, FUN=g, ...)
    val.sp <- gx.sp %*% weights
  }
  # The calculations in the above blocks apply the function g to every node
  # created then applies all the weights created
  
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

## Compare speed with parallel:

library(microbenchmark)
microbenchmark(sg.int(mySine, lower=rep(0, 3), upper=rep(9, 3), dimensions=3, parallel=FALSE),
               sg.int(mySine, lower=rep(0, 3), upper=rep(9, 3), dimensions=3, parallel=TRUE))
# For my computer, it was worse, but I only have one core

## Compare speed with adaptIntegrate:

library(cubature)
microbenchmark(sg.int(mySine, lower=rep(0, 2), upper=rep(9, 2), dimensions=2, parallel=FALSE),
               adaptIntegrate(mySine, rep(0, 2), rep(9, 2)))
# sg.int is about twice as fast



## Maximization:

func <- function(x, y){
  return(sin(x^2 / 2 - y^2 / 4) * cos(2*x - exp(y)))
}

y <- c(seq(from=-1, to=-0.01, by=0.01), seq(from=0.01, to=2.99, by=0.01))

a <- lapply(y, function(x) optimize(func, y=x, lower=c(-1, -1), upper=c(3, 3), maximum=TRUE))
# I couldn't quite get this one sorted

func <- function(x){
  return(sin(x[1]^2 / 2 - x[2]^2 / 4) * cos(2*x[1] - exp(x[2])))
}

optim(c(2, 1), func, method='L-BFGS-B', lower=c(-1, -1), upper=c(3, 3),
      control=list(fnscale=-1))
# but the output from this makes much more sense:
# The global max occurs around 2.03, 1.4