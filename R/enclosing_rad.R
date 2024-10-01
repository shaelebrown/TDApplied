#### COMPUTE enclosing RADIUS ####
#' Compute the enclosing radius for a dataset.
#'
#' The enclosing radius is the minimum radius beyond which no topological changes will occur.
#' 
#' @param X the input dataset, must either be a matrix or data frame.
#' @param distance_mat whether or not `X` is a distance matrix, default FALSE.
#' @return the numeric enclosing radius.
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @examples
#'
#' if(require("TDAstats"))
#' {
#'   # create a persistence diagram from a 2D Gaussian
#'   df = data.frame(x = rnorm(n = 20,mean = 0,sd = 1),y = rnorm(n = 20,mean = 0,sd = 1))
#' 
#'   # compute persistence diagram with calculate_homology from package TDAstats
#'   phom_TDAstats = TDAstats::calculate_homology(mat = df,dim = 0,threshold = 1)
#' 
#'   # convert to data frame
#'   phom_TDAstats_df = diagram_to_df(d = phom_TDAstats)
#' }
enclosing_radius <- function(X, distance_mat){
  
  # error check parameters
  if(is.null(distance_mat))
  {
    stop("distance_mat must not be NULL.")
  }
  if(length(distance_mat) > 1 | !inherits(distance_mat,"logical"))
  {
    stop("distance_mat must be a single logical (i.e. T or F).")
  }
  if(is.na(distance_mat) | is.nan(distance_mat) )
  {
    stop("distance_mat must not be NA/NAN.")
  }
  
  if(!inherits(X,"data.frame") & !inherits(X,"matrix"))
  {
    stop("X must either be a dataframe or a matrix.")
  }
  if(nrow(X) < 2 | ncol(X) < 1)
  {
    stop("X must have at least two rows and one column.")
  }
  if(length(which(stats::complete.cases(X) == F)) > 0)
  {
    stop("X must not contain any missing values.")
  }
  if(distance_mat == T & (ncol(X) != nrow(X) | !inherits(X,"matrix")))
  {
    stop("if distance_mat is TRUE then X must be a square matrix.")
  }
  if((inherits(X,"matrix") & !inherits(X[1,1],"numeric")) | (inherits(X,"data.frame") & length(which(unlist(lapply(X,is.numeric)))) < ncol(X)))
  {
    stop("X must have only numeric entries.")
  }
  
  # if X is not a distance matrix, compute distance mat
  if(!distance_mat)
  {
    dist_X <- dist(X)
    n <- nrow(X)
    return(min(sapply(1:(n - 1),FUN = function(X){
      
      lower_bound <- n*(X - 1) - X*(X - 1)/2 + 1
      upper_bound <- lower_bound + n - X
      if(X == n - 1)
      {
        upper_bound <- upper_bound - 1
      }
      return(max(dist_X[lower_bound:upper_bound]))
      
    })))
  }
  
  enc_rad <- min(apply(X, MARGIN = 1L, max))
  return(enc_rad)
  
}
