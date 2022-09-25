#### PERSISTENCE FISHER KERNEL ####
#' Calculate persistence Fisher kernel value between a pair of persistence diagrams.
#'
#' Returns the persistence Fisher kernel value between a pair of persistence diagrams
#' in a particular homological dimension, each of which is either the output from a \code{\link{diagram_to_df}} 
#' function call or from a TDA/TDAstats homology calculation like \code{\link[TDA]{ripsDiag}} or \code{\link[TDAstats]{calculate_homology}}.
#'
#' The persistence Fisher kernel is calculated from the Fisher information metric according to the formula
#' \eqn{k_{PF}(D_1,D_2) = exp(-t*d_{FIM}(D_1,D_2))}, resembling a radial basis kernel for standard
#' Euclidean spaces.
#'
#' @param D1 the first persistence diagram, either the output of a TDA/TDAstats calculation like \code{\link[TDA]{ripsDiag}}/\code{\link[TDAstats]{calculate_homology}}, or \code{\link{diagram_to_df}}.
#' @param D2 the second persistence diagram, either the output of a TDA/TDAstats calculation like \code{\link[TDA]{ripsDiag}}/\code{\link[TDAstats]{calculate_homology}}, or \code{\link{diagram_to_df}}.
#' @param dim the non-negative integer homological dimension in which the distance is to be computed, default 0.
#' @param sigma a positive number representing the bandwidth for the Fisher information metric, default 1.
#' @param t a positive number representing the scale for the persistence Fisher kernel, default 1.
#'
#' @return the numeric kernel value.
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @references
#' Le T, Yamada M (2018). "Persistence fisher kernel: a riemannian manifold kernel for persistence diagrams." \url{https://proceedings.neurips.cc/paper/2018/file/959ab9a0695c467e7caf75431a872e5c-Paper.pdf}.
#' 
#' Murphy, K. "Machine learning: a probabilistic perspective", MIT press (2012).
#' @examples
#'
#' # create two diagrams
#' D1 <- TDA::ripsDiag(TDA::circleUnif(n = 20,r = 1),
#'                     maxdimension = 1,maxscale = 2)
#' D2 <- TDA::ripsDiag(TDA::sphereUnif(n = 20,d = 2,r = 1),
#'                     maxdimension = 1,maxscale = 2)
#' 
#' # calculate the kernel value between D1 and D2 with sigma = 2, t = 2 in dimension 1
#' diagram_kernel(D1,D2,dim = 1,sigma = 2,t = 2)
#' # calculate the kernel value between D1 and D2 with sigma = 2, t = 2 in dimension 0
#' diagram_kernel(D1,D2,dim = 0,sigma = 2,t = 2)

diagram_kernel <- function(D1,D2,dim = 0,sigma = 1,t = 1){
  
  # check kernel-specific parameter, other inputs are checked in distance calculation
  check_param("t",t,non_negative = T,positive = F,numeric = T,finite = T,multiple = F)
  
  # return kernel calculation
  return(exp(-1*t*diagram_distance(D1 = D1,D2 = D2,dim = dim,distance = "fisher",sigma = sigma)))
  
}

#### GRAM MATRIX ####
#' Compute the gram matrix for a group of persistence diagrams.
#' 
#' Calculate the Gram matrix \eqn{K} for either a single list of persistence diagrams \eqn{(D_1,D_2,\dots,D_n)}, i.e. \eqn{K[i,j] = k_{PF}(D_i,D_j)}, 
#' or between two lists of persistence diagrams, \eqn{(D_1,D_2,\dots,D_n)} and \eqn{(D'_1,D'_2,\dots,D'_n)}, \eqn{K[i,j] = k_{PF}(D_i,D'_j)}, in parallel.
#' 
#' Gram matrices are used in downstream analyses, like in the `diagram_kkmeans`, `diagram_nearest_cluster`,`diagram_kpca`, 
#' `predict_diagram_kpca`, `predict_diagram_ksvm` and `independence_test` functions.
#'
#' @param diagrams a list of persistence diagrams, where each diagram is either the output of a TDA/TDAstats calculation like \code{\link[TDA]{ripsDiag}}/\code{\link[TDAstats]{calculate_homology}}, or \code{\link{diagram_to_df}}.
#' @param other_diagrams either NULL (default) or another list of persistence diagrams to compute a cross-Gram matrix.
#' @param dim the non-negative integer homological dimension in which the distance is to be computed, default 0.
#' @param sigma a positive number representing the bandwidth for the Fisher information metric, default 1.
#' @param t a positive number representing the scale for the kernel, default 1.
#' @param num_workers the number of cores used for parallel computation, default is one less than the number of cores on the machine.
#'
#' @return the numeric (cross) Gram matrix of class 'kernelMatrix'.
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @examples
#'
#' # create two diagrams
#' D1 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' D2 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                    dim = 0,threshold = 2)
#' g <- list(D1,D2)
#'
#' # calculate the Gram matrix in dimension 0 with sigma = 2, t = 2
#' G <- gram_matrix(diagrams = g,dim = 0,sigma = 2,t = 2,num_workers = 2)
#' 
#' # calculate cross-Gram matrix, which is the same as G
#' G_cross <- gram_matrix(diagrams = g,other_diagrams = g,dim = 0,sigma = 2,
#'                        t = 2,num_workers = 2)

gram_matrix <- function(diagrams,other_diagrams = NULL,dim = 0,sigma = 1,t = 1,num_workers = parallelly::availableCores(omit = 1)){
  
  # compute gram matrix from distance matrix
  K <- exp(-t*distance_matrix(diagrams = diagrams,other_diagrams = other_diagrams,dim = dim,distance = "fisher",sigma = sigma,num_workers = num_workers))
  
  # update class for interfacing with kernlab package
  class(K) <- "kernelMatrix"
  
  return(K)

}
