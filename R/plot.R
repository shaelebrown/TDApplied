#### PLOTTING PERSISTENCE DIAGRAMS ####
#' Plot persistence diagrams
#'
#' Plots a persistence diagram outputted from either a persistent homology calculation or from diagram_to_df, with
#' maximum homological dimension no more than 12 (otherwise the legend doesn't fit in the plot).
#' Each homological dimension has its own color and point type (with colors chosen to be clear and distinct from each other), 
#' and the main plot title can be altered via the `title` parameter.
#' 
#' @param D a persistence diagram, either outputted from either a persistent homology homology calculation like \code{\link[TDA]{ripsDiag}}/\code{\link[TDAstats]{calculate_homology}}/\code{\link{PyH}} or from \code{\link{diagram_to_df}}, with
#' maximum dimension at most 12.
#' @param title the character string plot title, default NULL.
#' @param max_radius the x and y limits of the plot are defined as `c(0,max_radius)`, and the default value of `max_radius` is the maximum death value in `D`.
#' @param legend a logical indicating whether to include a legend of feature dimensions, default TRUE.
#' @importFrom graphics legend abline
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @examples
#'
#' # create a sample diagram
#' diag <- TDAstats::calculate_homology(data.frame(x = rnorm(20),y = rnorm(20)))
#' 
#' # plot without title
#' plot_diagram(diag)
#' 
#' # plot with title
#' plot_diagram(diag,title = "Example diagram")
plot_diagram <- function(D,title = NULL,max_radius = NULL,legend = TRUE){
  
  # error check parameters
  check_diagram(d = D,ret = F)
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
  
  # convert diagram to df
  D <- diagram_to_df(D)
  
  # subset by max_radius
  D <- D[which(D$birth < max_radius & D$death <= max_radius),]
  
  # build plot
  pchs <- c(15:18,0:8)
  cols <- c("black","red","blue","darkgreen","violet","salmon","purple","orange","maroon","gold","chocolate","aquamarine","brown")
  dims <- unique(D[,1L])
  dims <- dims[order(dims)]
  plot(x = D[,2L],y = D[,3L],xlim = c(0,max_radius),ylim = c(0,max_radius),
       xlab = "Feature birth",ylab = "Feature death",col = cols[D[,1L] + 1],
       pch = pchs[D[,1L] + 1],main = ifelse(test = is.null(title),yes = "",no = title))
  if(legend == T)
  {
    graphics::legend("bottomright",title = "Dimensions:",legend = as.character(dims),col = cols[dims + 1],pch = pchs[dims + 1],inset = 0.01*max(D[,3L])) 
  }
  graphics::abline(a = 0,b = 1)
  
}

