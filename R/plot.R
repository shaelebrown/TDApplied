#### PLOTTING PERSISTENCE DIAGRAMS ####
#' Plot persistence diagrams
#'
#' Plots a persistence diagram outputted from either a persistent homology calculation or from diagram_to_df, with
#' maximum homological dimension no more than 12 (otherwise the legend doesn't fit in the plot).
#' Each homological dimension has its own color (the rcartocolor color-blind safe color palette) and point type, 
#' and the main plot title can be altered via the `title` parameter. Each feature is plotted with
#' a black point at its center in order to distinguish between overlapping features and easily compare
#' features to their persistence thresholds.
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
#' @importFrom graphics legend lines points
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @examples
#'
#' if(require("TDA") & require("TDAstats"))
#' {
#'   # create a sample diagram from the unit circle
#'   df <- TDA::circleUnif(n = 50)
#'   diag <- TDAstats::calculate_homology(df,threshold = 2)
#' 
#'   # plot without title
#'   plot_diagram(diag)
#' 
#'   # plot with title
#'   plot_diagram(diag,title = "Example diagram")
#' 
#'   # determine persistence thresholds
#'   thresholds <- bootstrap_persistence_thresholds(X = df,maxdim = 1,
#'   thresh = 2,num_samples = 3,
#'   num_workers = 2)
#' 
#'   # plot with bootstrap persistence thresholds
#'   plot_diagram(diag,title = "Example diagram with thresholds",thresholds = thresholds)
#' 
#'   #' # plot with personalized persistence thresholds
#'   plot_diagram(diag,title = "Example diagram with personalized thresholds",thresholds = c(0.5,1))
#' }

plot_diagram <- function(D,title = NULL,max_radius = NULL,legend = TRUE,thresholds = NULL){
  
  # error check parameters
  check_diagram(d = D,ret = F)
  
  # convert diagram to df
  D <- diagram_to_df(D)
  
  # error check title
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
    if(!inherits(title,"character"))
    {
      stop("title must be a character string.")
    }
  }
  
  # build plot
  pchs <- c(15:18,0:8)
  cols <- c("#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#88CCEE", 
            "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

  # if D non-empty, plot
  # otherwise plot empty diagram
  if(nrow(D) > 0)
  {
    if(max(D[,1L]) > 12)
    {
      stop("diagram must have maximum dimension no more than 12.")
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
    
    dims <- unique(D[,1L])
    dims <- dims[order(dims)]
    
    if(is.null(thresholds))
    {
      C <- cols[D[,1L] + 1]
    }else
    {
      C <- unlist(lapply(X = 1:nrow(D),FUN = function(X){
        
        return(ifelse(D[X,3L] - D[X,2L] > thresholds[[D[X,1L] + 1]],yes = cols[D[X,1L] + 1],no = "gray"))
        
      }))
    }
    if(is.null(thresholds))
    {
      plot(x = D[,2L],y = D[,3L],xlim = c(0,max_radius),ylim = c(0,max_radius),
           xlab = "Birth",ylab = "Death",col = C,
           pch = pchs[D[,1L] + 1],main = ifelse(test = is.null(title),yes = "",no = title))
      graphics::points(x = D[,2L],y = D[,3L],pch = pchs[[2]],col = "black",cex = 0.2)
    }else
    {
      gray_inds <- which(C == "gray")
      non_gray_inds <- which(C != "gray")
      plot(x = D[gray_inds,2L],y = D[gray_inds,3L],xlim = c(0,max_radius),ylim = c(0,max_radius),
           xlab = "Birth",ylab = "Death",col = C[gray_inds],
           pch = pchs[D[gray_inds,1L] + 1],main = ifelse(test = is.null(title),yes = "",no = title))
      graphics::points(x = D[non_gray_inds,2L],y = D[non_gray_inds,3L],col = C[non_gray_inds],pch = pchs[D[non_gray_inds,1L] + 1])
      graphics::points(x = D[,2L],y = D[,3L],pch = pchs[[2]],col = "black",cex = 0.2)
    }
    if(legend == T)
    {
      l <- c(expression(H[0]~' clusters'),expression(H[1]~' loops'),expression(H[2]~' voids'),expression(H[3]),expression(H[4]),expression(H[5]),expression(H[6]),expression(H[7]),expression(H[8]),expression(H[9]),expression(H[10]),expression(H[11]),expression(H[12]))[dims + 1]
      graphics::legend("bottomright",legend = l,col = cols[dims + 1],pch = pchs[dims + 1],inset = 0.01*max(D[,3L])) 
    }
    graphics::lines(x = c(0,max_radius),y = c(0,max_radius))
    if(!is.null(thresholds))
    {
      for(i in dims)
      {
        graphics::lines(x = c(0,max_radius - thresholds[[i + 1]]),y = c(thresholds[[i + 1]],max_radius),col = cols[i + 1],lty = "dashed")
      } 
    }
  }else
  {
    if(is.null(max_radius))
    {
      max_radius <- 1
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
    
    plot(x = numeric(),y = numeric(),xlim = c(0,max_radius),ylim = c(0,max_radius),
         xlab = "Feature birth",ylab = "Feature death",
         main = ifelse(test = is.null(title),yes = "",no = title))
    graphics::lines(x = c(0,max_radius),y = c(0,max_radius))
  }
  
}

