#### PLOTTING PERSISTENCE DIAGRAMS ####
#' Plot persistence diagrams
#'
#' Plots a persistence diagram outputted from either a persistent homology calculation or from diagram_to_df, with
#' maximum homological dimension no more than 3.
#' Each homological dimension has its own color and point type (with colors chosen to be clear and distinct from each other), 
#' and the main plot title can be altered via the `title` parameter.
#' 
#' @param D a persistence diagram, either outputted from either a persistent homology homology calculation like \code{\link[TDA]{ripsDiag}}/\code{\link[TDAstats]{calculate_homology}}/\code{\link{PyH}} or from \code{\link{diagram_to_df}}, with
#' maximum dimension at most 3.
#' @param title the character string plot title, default NULL.
#' @param xlim the x-limits of the plot as a length 2 numeric vector, default is c(0,max death value).
#' @param ylim the y-limits of the plot as a length 2 numeric vector, default is c(0,max death value).
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
plot_diagram <- function(D,title = NULL,xlim = NULL,ylim = NULL,legend = TRUE){
  
  # error check parameters
  check_diagram(d = D,ret = F)
  if(max(D[,1L]) > 3)
  {
    stop("diagram must have maximum dimension no more than 3.")
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
  if(is.null(xlim))
  {
    xlim <- c(0,max(D[,3L]))
  }else
  {
    if(!is.vector(xlim) | !is.numeric(xlim))
    {
      stop("xlim must be a numeric vector")
    }
    if(length(xlim) != 2)
    {
      stop("xlim must be a length 2 vector.")
    }
    if(xlim[[1]] >= xlim[[2]])
    {
      stop("lower limit of xlim must be smaller than its upper limit.")
    }
    if(is.infinite(xlim[[2]]))
    {
      stop("xlim must have a finite upper limit.")
    }
  }
  if(is.null(ylim))
  {
    ylim <- c(0,max(D[,3L]))
  }else
  {
    if(!is.vector(ylim) | !is.numeric(ylim))
    {
      stop("ylim must be a numeric vector")
    }
    if(length(ylim) != 2)
    {
      stop("ylim must be a length 2 vector.")
    }
    if(ylim[[1]] >= ylim[[2]])
    {
      stop("lower limit of ylim must be smaller than its upper limit.")
    }
    if(is.infinite(ylim[[2]]))
    {
      stop("ylim must have a finite upper limit.")
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
  
  # build plot
  pchs <- c(15:18)
  cols <- c("black","red","blue","darkgreen")
  dims <- unique(D[,1L])
  dims <- dims[order(dims)]
  plot(x = D[,2L],y = D[,3L],xlim = xlim,ylim = ylim,
       xlab = "Feature birth",ylab = "Feature death",col = cols[D[,1L] + 1],
       pch = pchs[D[,1L] + 1],main = ifelse(test = is.null(title),yes = "",no = title))
  if(legend == T)
  {
    graphics::legend("bottomright",title = "Dimensions:",legend = as.character(dims),col = cols[dims + 1],pch = pchs[dims + 1],inset = 0.01*max(D[,3L])) 
  }
  graphics::abline(a = 0,b = 1)
  
}

