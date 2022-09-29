#### PLOTTING PERSISTENCE DIAGRAMS ####
#' Plot persistence diagrams
#'
#' Plots a persistence diagram outputted from either a TDA/TDAstats homology calculation or from diagram_to_df, with
#' maximum homological dimension no more than 3.
#' Each homological dimension has its own color and point type (with colors chosen to be clear and distinct from each other), 
#' and the main plot title can be altered via the `title` parameter.
#' 
#' @param D a persistence diagram, either outputted from either a TDA/TDAstats homology calculation or from diagram_to_df, with
#' maximum dimension at most 3.
#' @param title the character string plot title, default NULL.
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
plot_diagram <- function(D,title = NULL){
  
  # error check parameters
  check_diagram(d = D)
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
  
  # convert diagram to df
  D <- diagram_to_df(D)
  
  # build plot
  pchs <- c(15:18)
  cols <- c("black","red","blue","darkgreen")
  dims <- unique(D[,1L])
  dims <- dims[order(dims)]
  plot(x = D[,2L],y = D[,3L],xlim = c(0,max(D[,3L])),ylim = c(0,max(D[,3L])),
       xlab = "Feature birth",ylab = "Feature death",col = cols[D[,1L] + 1],
       pch = pchs[D[,1L] + 1],main = ifelse(test = is.null(title),yes = "",no = title))
  abline(a = 0,b = 1)
  legend("bottomright",title = "Dimensions:",legend = as.character(dims),col = cols[dims + 1],pch = pchs[dims + 1],inset = 0.01*max(D[,3L]))
  
}