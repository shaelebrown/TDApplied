
#### ANALYZING MULTIPLE REPRESENTATIVES ####
#' Analyze the data point memberships of multiple representative (co)cycles.
#'
#' Multiple distance matrices with corresponding data points can contain the same topological features. 
#' Therefore we may wish to compare many representative (co)cycles across distance matrices to decide if their topological features are the same.
#' The `analyze_representatives` function returns a matrix of binary datapoint memberships in an input list of representatives across distance matrices.
#' Optionally this matrix can be plotted as a heatmap with columns as data points and rows (i.e. representatives) reordered by similarity, and the 
#' contributions (i.e. percentage membership) of each point in the representatives can also be returned. The heatmap has
#' dark red squares representing membership - location [i,j] is dark red if data point j is in representative i.
#' 
#' The clustering dendrogram can be used to determine if there are any similar groups of representatives (i.e.
#' shared topological features across datasets) and if so how many. The row labels of the heatmap are of the form
#' 'DX[Y]', meaning the Yth representative of diagram X, and the column labels are the data point numbers.
#' If diagrams are the output of the \code{\link{bootstrap_persistence_thresholds}}
#' function, then the subsetted_representatives (if present) will be analyzed. Therefore, a column label like 'DX[Y]' in the 
#' plotted heatmap would mean the Yth representative of diagram X. If certain representatives should be highlighted (by drawing a box around its row)
#' in the heatmap, a dataframe `boxed_reps` can be supplied with two integer columns - 'diagram' and 'rep'. For example, if we wish to draw a box for DX[Y] then we
#' add the row (diagram = X,rep = Y) to `boxed_reps`. If `d` is supplied then it will be used to cluster the representatives, based on the distances in `d`.
#'
#' @param diagrams a list of persistence diagrams, either the output of persistent homology calculations like ripsDiag/\code{\link[TDAstats]{calculate_homology}}/\code{\link{PyH}}, \code{\link{diagram_to_df}} or \code{\link{bootstrap_persistence_thresholds}}.
#' @param dim the integer homological dimension of representatives to consider.
#' @param num_points the integer number of data points in all the original datasets (from which the diagrams were calculated).
#' @param plot_heatmap a boolean representing if a heatmap of data point membership similarity of the representatives should be plotted, default `TRUE`. A dendrogram of hierarchical clustering is plotted, and rows (representatives) are sorted according to this clustering.
#' @param return_contributions a boolean indicating whether or not to return the membership contributions (i.e. percentages) of the data points (1:`num_points`) across all the representatives, default `FALSE`.
#' @param boxed_reps a data frame specifying specific rows of the output heatmap which should have a box drawn around them (for highlighting), default NULL. See the details section for more information.
#' @param lwd a positive number width for the lines of drawn boxes, if boxed_reps is not null.
#' @param d either NULL (default) or a "dist" object representing a distance matrix for the representatives, which must have the same number of rows and columns as cycles in the dimension `dim`.
#' @param title a character string title for the plotted heatmap, default NULL.
#' @param return_clust a boolean determining whether or not to return the result of the `stats::hclust` call when a heatmap is plotted, default `FALSE`.
#' @return either a matrix of data point contributions to the representatives, or a list with elements "memberships" (the matrix) and some combination of elements "contributions" (a vector of membership percentages for each data point across representatives) and "clust" (the results of `stats::hclust` on the membership matrix).
#' @importFrom stats heatmap order.dendrogram as.dendrogram hclust as.dist
#' @importFrom graphics rect
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#'   

analyze_representatives <- function(diagrams,dim,num_points,plot_heatmap = TRUE,return_contributions = FALSE,boxed_reps = NULL,d = NULL,lwd = NULL,title = NULL,return_clust = FALSE){
  
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
  if(length(plot_heatmap) > 1 | !inherits(plot_heatmap,"logical"))
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
  if(length(return_contributions) > 1 | !inherits(return_contributions,"logical"))
  {
    stop("return_contributions must be a single logical (i.e. T or F).")
  }
  if(is.na(return_contributions) | is.nan(return_contributions) )
  {
    stop("return_contributions must not be NA/NAN.")
  }
  if(plot_heatmap == T)
  {
    # error check title
    if(!is.null(title))
    {
      if(length(title) > 1)
      {
        stop("title must be a single character string.")
      }
      if(is.na(title) | is.nan(title))
      {
        stop("title must not be NA/NaN.")
      }
      if(!inherits(title,"character"))
      {
        stop("title must be a character string.")
      }
    }
  }
  if(is.null(return_clust))
  {
    stop("return_clust must not be NULL.")
  }
  if(length(return_clust) > 1 | !inherits(return_clust,"logical"))
  {
    stop("return_clust must be a single logical (i.e. T or F).")
  }
  if(is.na(return_clust) | is.nan(return_clust) )
  {
    stop("return_clust must not be NA/NAN.")
  }
  
  # check diagrams
  check_param("diagrams",diagrams,min_length = 2)
  diagrams <- lapply(X = diagrams,FUN = function(X){
    
    if("diagram" %in% names(X))
    {
      X$diag <- diagram_to_df(X)
    }
    
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
    if(!is.list(X$representatives[[1]]))
    {
      ripsDiag <- T
      # ripsDiag
      if(length(which(unlist(lapply(X = X$representatives,FUN = function(X){
        
        return(class(X))
        
      })) %in% c("matrix","array") == F)) > 0 | length(which(unlist(lapply(X = X$representatives,FUN = function(X){
        
        return(!all(X == round(X)))
        
      })))) > 0)
      {
        stop("Representatives must be matrices with integer entries - make sure that the diagrams were calculated from distance matrices.")
      }
    }else
    {
      # PyH
      ripsDiag <- F
      classes <- unlist(lapply(X = X$representatives,FUN = function(X){
        
        Y <- X
        if(length(Y) == 0)
        {
          return("array")
        }
        return(unlist(lapply(X = Y,FUN = function(X){
          
          if(!all(X == round(X)))
          {
            return("float")
          }
          return(class(X))
          
        })))
        
      }))
      if(all(classes %in% c("matrix","array")) == F)
      {
        stop("Representatives must be matrices with integer entries - make sure that the diagrams were calculated from distance matrices.")
      }
    }
    
    # make sure that all representatives do not have any data point greater than num_points
    if(ripsDiag)
    {
      if(suppressWarnings(max(unlist(lapply(X$representatives,FUN = max)))) > num_points)
      {
        stop("No representative should contain a data point which is greater than num_points.")
      }
    }else
    {
      if(suppressWarnings(max(unlist(lapply(X$representatives,FUN = function(X){
        
        Y <- X
        return(lapply(X = Y,FUN = function(X){
          
          return(max(X))
          
        }))
        
      })))) > num_points)
      {
        stop("No representative should contain a data point which is greater than num_points.")
      }
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
    check_diagram(X$diag,ret = F)
    
    return(list(diag = X$diag,representatives = X$representatives))
    
  })
  
  if(!is.null(boxed_reps))
  {
    if(inherits(boxed_reps,"data.frame") == F)
    {
      stop("If supplied, boxed_reps must be a data frame.")
    }
    if(ncol(boxed_reps) != 2 | all(colnames(boxed_reps) %in% c("diagram","rep")) == F)
    {
      stop("boxed_reps must have two columns - 'diagram' and 'rep'.")
    }
    if(nrow(boxed_reps) > 0)
    {
      if(length(which(boxed_reps$diagram != round(boxed_reps$diagram))) > 0 | length(which(boxed_reps$rep != round(boxed_reps$rep))) > 0)
      {
        stop("boxed_reps must contain only integer values.")
      }
      if(length(which(boxed_reps$diagram %in% 1:length(diagrams) == F)) > 0)
      {
        stop("boxed_reps contains a diagram entry outside of 1:length(diagrams).")
      }
    }else
    {
      boxed_reps <- NULL
    }
    
    if(!is.null(lwd))
    {
      check_param(param_name = "lwd",param = lwd,numeric = T,positive = T,multiple = F,finite = T) 
    }
    
  }
  
  cycle_inds <- lapply(X = diagrams,FUN = function(X){
    
    return(which(X$diag[,1] == dim))
    
  })
  
  cycles <- lapply(X = diagrams,FUN = function(X){
    
    lst <- X$representatives[which(X$diag[,1] == dim)]
    lst <- lapply(X = lst,FUN = function(X){
      
      return(unique(as.numeric(X)))
      
    })
    
  })
  
  if(!is.null(boxed_reps))
  {
    error_check <- lapply(X = 1:nrow(boxed_reps),FUN = function(X){
      
      if(boxed_reps[X,2L] %in% cycle_inds[[boxed_reps[X,1L]]] == F)
      {
        stop("For some row [i,j] of boxed_reps, the jth representative of diagram i did not have the input dimension 'dim'.")
      }
        
    }) 
  }
  
  num_cycles <- sum(unlist(lapply(X = cycles,FUN = length)))
  if(num_cycles == 0)
  {
    stop("No representatives were found in the desired dimension.")
  }
  
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
    if(!is.null(d))
    {
      if(!inherits(d,"dist"))
      {
        stop("If supplied, d must be of class 'dist'.")
      }
      if(length(d) != nrow(mat)*(nrow(mat) - 1)/2)
      {
        stop("d must have the same number of rows and columns as the number of cycles in the desired dimension.")
      }
    }else
    {
      d <- stats::dist(mat)
    }
    cluster <- stats::hclust(d = d)
    dendrogram  = stats::as.dendrogram(cluster)
    stats::heatmap(mat,scale = "none",xlab = "Data point",ylab = "Representative",main = title,Colv = NA,Rowv = dendrogram,add.expr = {
      
      if(!is.null(boxed_reps))
      {
        if(nrow(boxed_reps) > 0)
        {
          dend_order <- stats::order.dendrogram(dendrogram)
          box_top <- 0.5 + length(cycles)
          box_right <- 0.5 + num_points
          box_left <- 0.5
          box_bot <- 0.5
          box_height <- (box_top - box_bot)/nrow(mat)
          for(i in 1:nrow(boxed_reps))
          {
            ind <- which(rownames(mat) == paste0("D",boxed_reps[i,1L],"[",boxed_reps[i,2L],"]"))
            rownum <- which(dend_order == ind)
            if(!is.null(lwd))
            {
              graphics::rect(xleft = box_left,xright = box_right,ybottom = rownum - 0.5,ytop = rownum + 0.5,lwd = lwd,border = "red")
            }else
            {
              graphics::rect(xleft = box_left,xright = box_right,ybottom = rownum - 0.5,ytop = rownum + 0.5,border = "red")
            }
          }
        }
      }
      
    })
  }
  
  # compute contributions if desired
  if(return_contributions == T)
  {
    contributions <- apply(X = mat,MARGIN = 2L,mean)
    if(plot_heatmap == T & return_clust == T)
    {
      # return clustering object if desired
      return(list(mat = mat,contributions = contributions,clust = cluster))
    }else
    {
      return(list(mat = mat,contributions = contributions))
    }
  }else
  {
    if(plot_heatmap == T & return_clust == T)
    {
      # return clustering object if desired
      return(list(mat = mat,clust = cluster))
    }
  }
  
  return(mat)
  
}
