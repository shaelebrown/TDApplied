
#### ANALYZING MULTIPLE REPRESENTATIVES ####
#' Analyze the data point memberships of multiple representative (co)cycles.
#'
#' Multiple datasets with corresponding data points can contain the same topological features. 
#' Therefore we may wish to compare many representative (co)cycles across datasets to decide if their topological features are the same.
#' The `analyze_representatives` function returns a matrix of binary datapoint memberships in an input list of representatives across datasets.
#' Optionally this matrix can be plotted as a heatmap with rows (i.e. representatives) reordered by similarity, and the 
#' contributions (i.e. percentage membership) of each point in the representatives can also be returned.
#' 
#' The clustering dendrogram can be used to determine if there are any similar groups of representatives (i.e.
#' shared topological features across datasets) and if so how many.
#'
#' @param diagrams a list of persistence diagrams, either the output of persistent homology calculations like \code{\link[TDA]{ripsDiag}}/\code{\link[TDAstats]{calculate_homology}}/\code{\link{PyH}}, or \code{\link{diagram_to_df}}.
#' @param dim the integer homological dimension of representatives to consider.
#' @param num_points the integer number of data points in all the original datasets (from which the diagrams were calculated).
#' @param plot_heatmap a boolean representing if a heatmap of data point membership similarity of the representatives should be plotted, default `TRUE`. A dendrogram of hierarchical clustering is plotted, and rows (representatives) are sorted according to this clustering.
#' @param return_contributions a boolean indicating whether or not to return the membership contributions (i.e. percentages) of the data points (1:`num_points`) across all the representatives, default `FALSE`.
#' @return either a matrix of data point contributions to the representatives, or a list with elements "memberships" (the matrix) and "contributions" (a vector of membership percentages for each data point across representatives).
#' @export
#' @importFrom methods is
#' @importFrom stats heatmap
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @examples
#'
#' if(require("TDA"))
#' {
#'   # sample the unit circle
#'   df <- TDA::circleUnif(n = 50)
#'                         
#'   # create 10 copies with added Gaussian noise and
#'   # calculate their diagrams
#'   circs <- lapply(X = 1:10,FUN = function(X){
#'      df <- circ
#'      df$x <- df$x + rnorm(n = 50,sd = 0.05)
#'      df$y <- df$y + rnorm(n = 50,sd = 0.05)
#'      diag <- bootstrap_persistence_thresholds(X = as.matrix(dist(df)),
#'                                               FUN = "ripsDiag",maxdim = 1,
#'                                               thresh = 2,distance_mat = T,
#'                                               calculate_representatives = T,
#'                                               return_subsetted = T)
#'      return(diag)
#'   
#'    })
#'    
#'    # analyze loop representatives across the diagrams
#'    # num_points is 50 because each underlying dataset had
#'    # 50 (corresponding) points
#'    analyze_representatives(diagrams = circs,dim = 1,
#'                            num_points = 50)
#'   
#'  }

analyze_representatives <- function(diagrams,dim,num_points,plot_heatmap,return_contributions){
  
  # error check parameters
  check_param(param_name = "dim",dim,numeric = T,at_least_one = T,multiple = F,infinite = F)
  check_param(param_name = "num_points",num_points,numeric = T,at_least_one = T,multiple = F,infinite = F)
  if(num_points < 4)
  {
    stop("num_points must be at least 4.")
  }
  if(is.null(plot_heatmap))
  {
    stop("plot_heatmap must not be NULL.")
  }
  if(length(plot_heatmap) > 1 | !methods::is(plot_heatmap,"logical"))
  {
    stop("plot_heatmap must be a single logical (i.e. T or F).")
  }
  if(is.na(plot_heatmap) | is.nan(plot_heatmap) )
  {
    stop("plot_heatmap must not be NA/NAN.")
  }
  if(is.null(return_contributions))
  {
    stop("return_contributions must not be NULL.")
  }
  if(length(return_contributions) > 1 | !methods::is(return_contributions,"logical"))
  {
    stop("return_contributions must be a single logical (i.e. T or F).")
  }
  if(is.na(return_contributions) | is.nan(return_contributions) )
  {
    stop("return_contributions must not be NA/NAN.")
  }
  
  # check diagrams
  check_param("diagrams",diagrams,min_length = 2)
  diagrams <- lapply(X = diagrams,FUN = function(X){
    
    # make sure there are either cycles or representatives, and rename
    if("cycleLocation" %in% names(X) == F & "representatives" %in% names(X) == F)
    {
      stop("Representative cycles must be stored in the diagrams, either calculated with ripsDiag or PyH.")
    }
    if("cycleLocation" %in% names(X))
    {
      X$representatives <- X$cycleLocation
    }
    # make sure that the representatives are lists of integer matrices
    if(!is.list(X$representatives))
    {
      stop("Representatives must be a list.")
    }
    if(length(which(unlist(lapply(X = X$representatives,FUN = function(X){
      
      return(class(X))
      
    })) %in% c("matrix","array") == F)) > 0 | length(which(unlist(lapply(X = X$representatives,FUN = function(X){
      
      return(!all(X == round(X)))
      
    })))) > 0)
    {
      stop("Representatives must be matrices with integer entries - make sure that the diagrams were calculated from distance matrices.")
    }
    
    # make sure that all representatives do not have any data point greater than num_points
    if(max(X) > num_points)
    {
      stop("No representative should contain a data point which is greater than num_points.")
    }
    
    # check persistence diagrams and rename
    if("diagram" %in% names(X) == F & "diag" %in% names(X) == F & "subsetted_diag" %in% names(X) == F)
    {
      stop("Persistence diagrams must be supplied.")
    }
    if("subsetted_diag" %in% names(X))
    {
      X$diag <- X$subsetted_diag
      X$representatives <- X$subsetted_representatives
    }
    if("diagram" %in% names(X))
    {
      X$diag <- diagram_to_df(X)
    }
    check_diagram(X$diag,ret = F)
    
  })
  
  cycle_inds <- lapply(X = diagrams,FUN = function(X){
    
    return(which(X$diag[,1] == dim))
    
  })
  
  cycles <- lapply(X = diagrams,FUN = function(X){
    
    lst <- X$representatives[which(X$diag[,1] == dim)]
    lst <- lapply(X = lst,FUN = function(X){
      
      return(unique(as.numeric(X)))
      
    })
    
  })
  
  mat <- do.call(rbind,lapply(X = cycles,FUN = function(X){
    
    Y <- X
    rows <- do.call(rbind,lapply(X = Y,FUN = function(X){
      
      v <- rep(0,num_points)
      v[X] <- 1
      return(v)
      
    }))
    
  }))
  rownames(mat) <- unlist(lapply(X = 1:length(diagrams),FUN = function(X){
    
    if(length(cycle_inds[[X]]) == 0)
    {
      return(c())
    }
    return(paste("D",X,"[",cycle_inds[[X]],"]",sep = ""))
    
  }))
  
  if(plot_heatmap)
  {
    stats::heatmap(mat,scale = "none",Colv = NA) 
  }
  
  if(return_contributions == T)
  {
    return(list(mat = mat,contributions = apply(X = mat,MARGIN = 2L,mean)))
  }
  
  return(mat)
  
}


