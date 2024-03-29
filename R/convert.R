#### CONVERT PERSISTENCE DIAGRAMS INTO DATA FRAMES####
#' Convert a TDA/TDAstats persistence diagram to a data frame.
#'
#' The output of homology calculations from the R packages TDA
#' and TDAstats are not dataframes. This function converts these 
#' outputs into a data frame either for further usage in this package or
#' for personalized analyses.
#' 
#' If a diagram is constructed using a TDA function like ripsDiag
#' with the `location` parameter set to true then the return value will ignore the location information.
#'
#' @param d the output of a TDA/TDAstats homology calculation, like ripsDiag or \code{\link[TDAstats]{calculate_homology}}.
#' @return a 3-column data frame, with each row representing a topological feature. The first column is the feature dimension (a non-negative integer), the second column is the birth radius of the feature and the third column is the death radius.
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

diagram_to_df <- function(d){

  # function to convert d to a data frame with standardized column names
  # d is a diagram from library TDA or TDAstats
  
  # preliminary check, mostly for internal methods
  if(inherits(d,"data.frame"))
  {
    return(d)
  }
  
  if((is.list(d) && ((length(d) == 1 && all(names(d) %in% "diagram") && (inherits(d$diagram,"diagram")) || inherits(d$diagram,"data.frame")) || ((length(d) == 4 && all(names(d) %in% c("diagram","birthLocation","deathLocation","cycleLocation")) && inherits(d$diagram,"diagram"))))) == F && (inherits(d,"matrix") && inherits(d,"array") & all(colnames(d) %in% c("dimension","birth","death"))) == F)
  {
    stop("Diagrams must either be the output of a TDA/TDAstats/PyH computation.")
  }
  
  if(inherits(d,"matrix") & inherits(d,"array"))
  {
    # diagram was the output of a TDAstats calculation
    return(as.data.frame(d))
  }
  
  if("diagram" %in% names(d))
  {
    if(inherits(d$diagram,"data.frame"))
    {
      # diagram was the output of a PyH calculation, with representatives
      return(d$diagram)
    }
  }

  # else d was the output of a TDA calculation
  d <- d[[1]]
  class(d) <- "matrix"
  d <- as.data.frame(d)
  colnames(d) <- c("dimension","birth","death")

  return(d)

}
