#### CONVERT TDA PERSISTENCE DIAGRAMS ####
#' Convert TDA persistence diagrams to data frames
#'
#' The output of homology calculations in package TDA
#' are not immediately able to be used for distance calculations.
#' This function converts the output of homology calculations in TDA
#' into a data frame either for further usage in this package or
#' for personalized analyses.
#'
#' The `d` parameter should be the output of a TDA homology calculation.
#'
#' @param d the output of a TDA homology calculation, like ripsDiag.
#' @return a 3-column data frame, with each row representing a topological feature.
#' @export
#' @examples
#'
#' # create a persistence diagram from a 2D Gaussian using TDA
#' df = data.frame(x = rnorm(n = 100,mean = 0,sd = 1),y = rnorm(n = 100,mean = 0,sd = 1))
#'
#' # compute persistence diagram with ripsDiag
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
