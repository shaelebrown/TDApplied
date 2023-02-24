
# functions to calculate Fisher information distance matrices and Gram matrices
# in parallel with a fast approximation

# these matrices can then be input into TDApplied functions directly

parallel_approx_distance_matrix <- function(diagrams,other_diagrams = NULL,dim = 0,sigma = 1,rho = 1e-3,num_workers = parallelly::availableCores(omit = 1)){
  
  # create cluster
  cl <- parallel::makeCluster(num_workers)
  doParallel::registerDoParallel(cl)
  
  # export TDApplied, and input parameters, to all cluster workers
  parallel::clusterEvalQ(cl,library("TDApplied"))
  parallel::clusterExport(cl,varlist = c("diagrams","dim","sigma","rho"),envir = environment())
  
  # calculate distances in parallel
  # clusters are closed
  tryCatch(expr = {
    
    if(is.null(other_diagrams))
    {
      # not cross distance matrix, only need to compute the upper diagonal
      # since the matrix is symmetric
      d <- matrix(data = 0,nrow = length(diagrams),ncol = length(diagrams))
      d_off_diag <- foreach::`%dopar%`(obj = foreach::foreach(r = iterators::iter(which(upper.tri(d),arr.ind = T),by = 'row'),.combine = c),ex = {TDApplied::diagram_distance(D1 = diagrams[[r[[1]]]],D2 = diagrams[[r[[2]]]],dim = dim,distance = "fisher",sigma = sigma,rho = rho)})
      d[upper.tri(d)] <- d_off_diag
      d[which(upper.tri(d),arr.ind = T)[,c("col","row")]] <- d_off_diag
      diag(d) <- rep(0,nrow(d))
    }else
    {
      # cross distance matrix, need to compute all entries
      d <- foreach::`%dopar%`(foreach::`%:%`(foreach::foreach(r = 1:length(other_diagrams),.combine = cbind),foreach::foreach(X = 1:length(diagrams),.combine = c)),ex = {TDApplied::diagram_distance(D1 = other_diagrams[[r]],D2 = diagrams[[X]],dim = dim,distance = "fisher",sigma = sigma,rho = rho)})
    }
    
  }, warning = function(w){warning(w)},
  error = function(e){stop(e)},
  finally = {
    # close cluster
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
    
  })
  
  return(d)
  
}

parallel_approx_gram_matrix <- function(diagrams,other_diagrams = NULL,dim = 0,sigma = 1,t = 1,rho = 1e-3,num_workers = parallelly::availableCores(omit = 1)){
  
  # compute gram matrix from distance matrix
  K <- exp(-t*parallel_approx_distance_matrix(diagrams = diagrams,other_diagrams = other_diagrams,dim = dim,sigma = sigma,rho = rho,num_workers = num_workers))
  
  # update class for interfacing with kernlab package
  class(K) <- "kernelMatrix"
  
  return(K)
  
}
