#### CONVERT TDA PERSISTENCE DIAGRAMS ####
#' Convert a TDA persistence diagram to a data frame.
#'
#' The output of homology calculations from the R package TDA
#' are not immediately able to be used for distance calculations.
#' This function converts the output of homology calculations in TDA
#' into a data frame either for further usage in this package or
#' for personalized analyses.
#' 
#' If a diagram is constructed using a TDA function like \code{\link[TDA]{ripsDiag}}
#' with the `location` parameter set to true then the return value will ignore the location information.
#'
#' @param d the output of a TDA homology calculation, like \code{\link[TDA]{ripsDiag}}.
#' @return a 3-column data frame, with each row representing a topological feature. The first column is the feature dimension (a non-negative integer), the second column is the birth radius of the feature and the third column is the death radius.
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @examples
#'
#' # create a persistence diagram from a 2D Gaussian
#' df = data.frame(x = rnorm(n = 100,mean = 0,sd = 1),y = rnorm(n = 100,mean = 0,sd = 1))
#'
#' # compute persistence diagram with ripsDiag from package TDA
#' phom = TDA::ripsDiag(X = df,maxdimension = 1,maxscale = 1)
#'
#' # convert to data frame
#' phom_df = diagram_to_df(d = phom)

diagram_to_df <- function(d){

  # function to convert d to a data frame with standardized column names
  # d is a diagram from library TDA

  d <- d[[1]]
  class(d) <- "matrix"
  d <- as.data.frame(d)
  colnames(d) <- c("dimension","birth","death")

  return(d)

}
