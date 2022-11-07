#### PLOTTING PERSISTENCE DIAGRAMS ####
#' Plot persistence diagrams
#'
#' Plots a persistence diagram outputted from either a persistent homology calculation or from diagram_to_df, with
#' maximum homological dimension no more than 12 (otherwise the legend doesn't fit in the plot).
#' Each homological dimension has its own color and point type (with colors chosen to be clear and distinct from each other), 
#' and the main plot title can be altered via the `title` parameter.
#' 
#' The `thresholds` parameter, if not NULL, can either be a user-defined numeric vector, with
#' one entry (persistence threshold) for each dimension in `D`, or the output of
#' \code{\link{bootstrap_persistence_thresholds}}. Points whose persistence are greater than or equal to their dimension's
#' threshold will be plotted in their dimension's color, and in gray otherwise.
#' 
#' @param D a persistence diagram, either outputted from either a persistent homology homology calculation like \code{\link[TDA]{ripsDiag}}/\code{\link[TDAstats]{calculate_homology}}/\code{\link{PyH}} or from \code{\link{diagram_to_df}}, with
#' maximum dimension at most 12.
#' @param title the character string plot title, default NULL.
#' @param max_radius the x and y limits of the plot are defined as `c(0,max_radius)`, and the default value of `max_radius` is the maximum death value in `D`.
#' @param legend a logical indicating whether to include a legend of feature dimensions, default TRUE.
#' @param thresholds either a numeric vector with one persistence threshold for each dimension in `D` or the output of a \code{\link{bootstrap_persistence_thresholds}} function call, default NULL.
#' @importFrom graphics legend abline
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @examples
#'
#' # create a sample diagram from the unit circle
#' df <- TDA::circleUnif(n = 50)
#' diag <- TDAstats::calculate_homology(df,threshold = 2)
#' 
#' # plot without title
#' plot_diagram(diag)
#' 
#' # plot with title
#' plot_diagram(diag,title = "Example diagram")
#' 
#' # determine persistence thresholds
#' thresholds <- bootstrap_persistence_thresholds(X = df,maxdim = 1,
#' thresh = 2,num_samples = 3,
#' num_workers = 2)
#' 
#' # plot with bootstrap persistence thresholds
#' plot_diagram(diag,title = "Example diagram with thresholds",thresholds = thresholds)
#' 
#' #' # plot with personalized persistence thresholds
#' plot_diagram(diag,title = "Example diagram with personalized thresholds",thresholds = c(0.5,1))
plot_diagram <- function(D,title = NULL,max_radius = NULL,legend = TRUE,thresholds = NULL){
  
  # error check parameters
  check_diagram(d = D,ret = F)
  
  # convert diagram to df
  D <- diagram_to_df(D)
  
  if(max(D[,1L]) > 12)
  {
    stop("diagram must have maximum dimension no more than 12.")
  }
  if(!is.null(title))
  {
    if(is.na(title) | is.nan(title))
    {
      stop("title must not be NA/NaN.")
    }
    if(length(title) > 1)
    {
      stop("title must be a single character string.")
    }
    if(!methods::is(title,"character"))
    {
      stop("title must be a character string.")
    }
  }
  if(is.null(max_radius))
  {
    max_radius <- max(D[,3L])
    if(is.infinite(max_radius))
    {
      stop("when D contains Inf death values either remove those points, set their death value to be finite, or set a finite max_radius value.")
    }
  }else
  {
    if(!is.vector(max_radius) | !is.numeric(max_radius))
    {
      stop("max_radius must be numeric.")
    }
    if(length(max_radius) != 1)
    {
      stop("max_radius must be a single number.")
    }
    if(is.infinite(max_radius))
    {
      stop("max_radius must be finite.")
    }
    if(max_radius <= 0)
    {
      stop("max_radius must be positive.")
    }
  }
  if(is.null(legend))
  {
    stop("legend must not be NULL.")
  }
  if(length(legend)>1 | !is.logical(legend))
  {
    stop("legend must be a single logical value.")
  }
  
  # error check thresholds
  if(!is.null(thresholds))
  {
    # first see if thresholds was the output of the bootstrap
    # function with multiple return arguments
    bootstrap_output <- F
    if(is.list(thresholds))
    {
      if("thresholds" %in% names(thresholds))
      {
        # if yes, then subset for just the thresholds
        thresholds = thresholds$thresholds
        bootstrap_output <- T
      }else
      {
        stop("thresholds must contain a list element called \'thresholds\'.")
      }
    }
    if(!is.vector(thresholds) | !is.numeric(thresholds))
    {
      if(bootstrap_output)
      {
        stop("the \'thresholds\' element of thresholds should be a numeric vector.")
      }else
      {
        stop("thresholds should be a numeric vector.")
      }
    }
    if(length(thresholds) != length(unique(D[,1])))
    {
      if(bootstrap_output)
      {
        stop("the \'thresholds\' element of thresholds must have one element for each dimension in D.")
      }else
      {
        stop("thresholds must have one element for each dimension in D.")
      }
    }
    if(length(which(thresholds %in% c(NA,NaN,Inf))) > 0)
    {
      if(bootstrap_output)
      {
        stop("the \'thresholds\' element of thresholds must not have any NA, NaN or Inf values.") 
      }else
      {
        stop("thresholds must not have any NA, NaN or Inf values.")
      }
    }
  }
  
  # subset by max_radius
  D <- D[which(D$birth < max_radius & D$death <= max_radius),]
  
  # build plot
  pchs <- c(15:18,0:8)
  cols <- c("black","red","blue","darkgreen","violet","salmon","purple","orange","maroon","gold","chocolate","aquamarine","brown")
  dims <- unique(D[,1L])
  dims <- dims[order(dims)]
  if(is.null(thresholds))
  {
    C <- cols[D[,1L] + 1]
  }else
  {
    C <- unlist(lapply(X = 1:nrow(D),FUN = function(X){
      
      return(ifelse(D[X,3L] - D[X,2L] >= thresholds[[D[X,1L] + 1]],yes = cols[D[X,1L] + 1],no = "gray"))
      
    }))
  }
  plot(x = D[,2L],y = D[,3L],xlim = c(0,max_radius),ylim = c(0,max_radius),
       xlab = "Feature birth",ylab = "Feature death",col = C,
       pch = pchs[D[,1L] + 1],main = ifelse(test = is.null(title),yes = "",no = title))
  if(legend == T)
  {
    graphics::legend("bottomright",title = "Dimensions:",legend = as.character(dims),col = cols[dims + 1],pch = pchs[dims + 1],inset = 0.01*max(D[,3L])) 
  }
  graphics::abline(a = 0,b = 1)
  if(!is.null(thresholds))
  {
    for(i in dims)
    {
      graphics::abline(b = 1,a = thresholds[[i + 1]],col = cols[i + 1],lty = "dashed")
    } 
  }
}

