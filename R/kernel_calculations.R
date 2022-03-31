#### PERSISTENCE FISHER KERNEL ####
#' Calculate persistence Fisher kernel between a pair of persistence diagrams
#'
#' Returns the persistence Fisher kenel between a pair of persistence diagrams
#' in a particular homological dimension.
#'
#' The `D1` and `D2` parameters are the two persistence diagrams.
#' The `dim` parameter should be a positive finite integer.
#' The `sigma` parameter is the positive bandwith for the Fisher information metric, and
#' `t` is the scale parameter for persistence Fisher kernel.
#'
#' @param D1 the first persistence diagram, either outputted from TDA or from diagram_to_df function.
#' @param D2 the second persistence diagram, either outputted from TDA or from diagram_to_df function.
#' @param dim the homological dimension in which the distance is to be computed.
#' @param sigma a positive number representing the bandwith for the Fisher information metric, default 1.
#' @param t a positive number representing the scale for the kernel, default 1.
#'
#' @return the kernel value
#' @export
#' @examples
#'
#' # create two diagrams with package TDA based on 2D Gaussians
#' diag1 <- TDA::ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),
#' y = rnorm(100,mean = 0,sd = 1)),
#' maxscale = 1,
#' maxdimension = 1)
#' diag2 <- TDA::ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),
#' y = rnorm(100,mean = 0,sd = 1)),
#' maxscale = 1,
#' maxdimension = 1)
#'
#' # calculate their kernel value in dimension 1 with sigma = 2, t = 2
#' k <- diagram_kernel(D1 = diag1,D2 = diag2,dim = 1,sigma = 2,t = 2)

diagram_kernel <- function(D1,D2,dim = 0,sigma = 1,t = 1){
  
  # check kernel-specific parameters, other inputs are checked in distance calculation
  if(is.null(t))
  {
    stop("t must not be NULL.")
  }
  if(!is.numeric(t) | length(t) > 1 | is.na(t) | is.nan(t) | t <= 0)
  {
    stop("t must be a positive number.")
  }
  
  # return kernel calculation
  return(exp(-1*t*diagram_distance(D1 = D1,D2 = D2,dim = dim,distance = "fisher",sigma = sigma)))
  
}

#### GRAM MATRIX ####
#' Calculate the Gram matrix K for a list of persistence diagrams, i.e. K[i,j] = k(d[i],d[j]), in parallel.
#'
#' Returns the Gram matrix for a list of persistence diagrams.
#'
#' `diagrams` is the list of persistence diagrams.
#' The `dim` parameter should be a positive finite integer.
#' The `sigma` parameter is the positive bandwith for the Fisher information metric, and
#' `t` is the positive scale parameter for the persistence Fisher kernel.
#'
#' @param diagrams the list of persistence diagrams, either the output from TDA calculations or the diagram_to_df function.
#' @param dim the homological dimension in which the distance is to be computed.
#' @param sigma a positive number representing the bandwith for the Fisher information metric, default 1.
#' @param t a positive number representing the scale for the kernel, default 1.
#'
#' @return the Gram matrix of class 'kernelMatrix' for downstream analyses with the kernlab package.
#' @export
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ
#' @importFrom parallelly availableCores
#' @importFrom doParallel registerDoParallel
#' @importFrom iterators iter
#' @examples
#'
#' # create ten diagrams with package TDA based on 2D Gaussians
#' g <- lapply(X = 1:10,FUN = function(X){
#'
#' diag <- TDA::ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),
#' y = rnorm(100,mean = 0,sd = 1)),
#' maxscale = 1,
#' maxdimension = 1)
#' df <- diagram_to_df(d = diag)
#' return(df)
#'
#' })
#'
#' # calculate their Gram matrix in dimension 1 with sigma = 2, t = 2
#' G <- gram_matrix(diagrams = g,dim = 1,sigma = 2,t = 2)

gram_matrix <- function(diagrams,dim = 0,sigma = 1,t = 1){
  
  # set internal variables to NULL to avoid build issues
  r <- NULL
  
  # error check diagrams argument
  if(is.null(diagrams))
  {
    stop("diagrams must be a list of persistence diagrams.")
  }
  if(!is.list(diagrams) | length(diagrams) < 2)
  {
    stop("diagrams must be a list of persistence diagrams of length at least 2.")
  }
  diagrams <- all_diagrams(diagram_groups = list(diagrams),inference = "independence")[[1]]
  
  # compute Gram matrix in parallel
  m = length(diagrams)
  num_workers <- parallelly::availableCores(omit = 1)
  cl <- parallel::makeCluster(num_workers)
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl,c("diagram_distance","diagram_kernel"))
  force(diagrams) # required for parallel computation in this environment
  force(check_diagram)
  
  K <- matrix(data = 1,nrow = m,ncol = m)
  k <- foreach::`%dopar%`(obj = foreach::foreach(r = iterators::iter(which(upper.tri(K),arr.ind = T),by = 'row'),.combine = c),ex = {
    
    return(diagram_kernel(D1 = diagrams[[r[[1]]]],D2 = diagrams[[r[[2]]]],dim = dim,sigma = sigma,t = t))
    
  })
  K[upper.tri(K)] <- k
  K[which(upper.tri(K),arr.ind = T)[,c("col","row")]] <- k
  
  parallel::stopCluster(cl)
  
  # update class for interfacing with kernlab package
  class(K) <- "kernelMatrix"
  
  return(K)

}



