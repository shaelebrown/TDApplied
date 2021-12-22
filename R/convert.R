#### CONVERT TDA PERSISTENCE DIAGRAMS ####
#' Convert TDA persistence diagrams to data frames
#'
#' The output of homology calculations in package TDA
#' are not immediately able to be used for distance calculations.
#' These functions convert the output of homology calculations in TDA
#' into a data frame for further usage.
#'
#' The `d` parameter is the persistence diagram to be converted.
#'
#' @param d output of a homology calculation, like ripsDiag
#' @return 3-column data frame, with each row representing a topological feature
#' @export
#' @examples
#'
#' # create a persistence diagram from a 2D Gaussian using TDA
#' df = data.frame(x = rnorm(n = 100,mean = 0,sd = 1),y = rnorm(n = 100,mean = 0,sd = 1))
#' # compute persistence barcode with ripsDiag
#' phom = ripsDiag(X = df,maxdimension = 1,maxscale = 1)
#' # convert to data frame
#' phom_df = TDA_diagram_to_df(d = phom)
#'
#' # create a persistence diagram from a 2D Gaussian using TDAStats
#' df = data.frame(x = rnorm(n = 100,mean = 0,sd = 1),y = rnorm(n = 100,mean = 0,sd = 1))
#'
#' # compute persistence barcode with calculate_homology
#' phom = calculate_homology(mat = df,maxdimension = 1,maxscale = 1)
#'
#' # convert to data frame
#' phom_df = TDAStats_diagram_to_df(d = phom,dim = 1,threshold = 1)

TDA_diagram_to_df <- function(d){

  # function to convert d to a data frame with standardized column names
  # d is a diagram from library TDA

  d <- d[[1]]
  class(d) <- "matrix"
  d <- as.data.frame(d)
  colnames(d) <- c("dimension","birth","death")

  return(d)

}

#### CONVERT TDAStats PERSISTENCE DIAGRAMS ####
#' Convert TDAStats persistence diagrams to data frames
#'
#' The output of homology calculations in package TDAStats
#' are not immediately able to be used for distance calculations.
#' These functions convert the output of homology calculations in TDAStats
#' into a data frame for further usage.
#'
#' The `d` parameter is the persistence diagram to be converted.
#'
#' @param d output of a calculate_homology calculation
#' @return 3-column data frame, with each row representing a topological feature
#' @export
#' @examples
#'
#' # create a persistence diagram from a 2D Gaussian using TDAStats
#' df = data.frame(x = rnorm(n = 100,mean = 0,sd = 1),y = rnorm(n = 100,mean = 0,sd = 1))
#'
#' # compute persistence barcode with calculate_homology
#' phom = calculate_homology(mat = df,maxdimension = 1,maxscale = 1)
#'
#' # convert to data frame
#' phom_df = TDAStats_diagram_to_df(d = phom,dim = 1,threshold = 1)

TDAStats_diagram_to_df <- function(d){

  # function to convert d to a data frame with standardized column names
  # d is a diagram from library TDAStats

  return(as.data.frame(d))

}
