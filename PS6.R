sg.int <- function(g, ..., lower, upper, dimensions){
  require("SparseGrid")
  require('parallel')
  
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
  
  clus <- makeCluster(max(detectCores() - 1, 1))
  gx.sp <- parApply(cl=clus, X=nodes, MARGIN=1, FUN=g, ...)
  val.sp <- gx.sp %*% weights
  stopCluster(clus)
  
  return(val.sp)
}