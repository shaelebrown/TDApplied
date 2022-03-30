#### PERSISTENCE FISHER KERNEL ####
#' Calculate persistence Fisher kernel between a pair of persistence diagrams
#'
#' Returns a function which calculates the persistence Fisher kenel between a pair of persistce diagrams, stored as
#' data frames (as the output from diagram_to_df) in a particular homological dimension.
#'
#' The `dim` parameter should be a positive finite integer.
#' The `sigma` parameter is the positive bandwith for the Fisher information metric, and
#' `t` is the scale parameter for persistence Fisher kernel.
#'
#' @param dim the homological dimension in which the distance is to be computed.
#' @param sigma a positive number representing the bandwith for the Fisher information metric, default 1.
#' @param t a positive number representing the scale for the kernel, default 1.
#'
#' @return the kernel function.
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
#' fisher_kernel <- diagram_kernel(dim = 1,sigma = 2,t = 2)
#' k <- fisher_kernel(D1 = diag1,D2 = diag2)

diagram_kernel <- function(dim = 0,sigma = 1,t = 1){
  
  # check kernel-specific parameters, other inputs are checked in distance calculation
  if(is.null(t))
  {
    stop("t must not be NULL.")
  }
  if(!is.numeric(t) | length(t) > 1 | is.na(t) | is.nan(t) | t <= 0)
  {
    stop("t must be a positive number.")
  }
  
  # return kernel function
  f <- function(D1,D2){
    
    return(exp(-1*t*diagram_distance(D1 = D1,D2 = D2,dim = dim,distance = "fisher",sigma = sigma)))
    
  }
  class(f) <- "kernel"
  
  return(f)
  
}





