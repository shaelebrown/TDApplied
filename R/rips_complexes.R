#### COMPUTE STATIC SIMPLICIAL COMPLEXES WITHIN RIPS FILTRATION ####
#' Compute Rips graphs of a dataset at particular epsilon radius values.
#'
#' Persistence diagrams computed from Rips-Vietoris filtrations contain information about 
#' distance radius scales at which topological features of a dataset exist, but the features
#' can be challenging to visualize, analyze and interpret. In order to help solve this problem the `rips_graphs`
#' function computes the 1-skeleton (i.e. graph) of Rips complexes at particular radii.
#' 
#' This function may be used in conjunction with the \link{igraph} package to visualize the graphs (see \code{\link{plot_rips_graph}}).
#'
#' @param X either a point cloud data frame/matrix, or a distance matrix.
#' @param distance_mat a boolean representing if the input `X` is a distance matrix, default value is `FALSE`.
#' @param eps a numeric vector of the positive scales at which to compute the Rips-Vietoris complexes, i.e. all edges at most the specified values.
#' @param return_clusters a boolean determining if the connected components (i.e. data clusters) of the complex should be explicitly returned, default is `TRUE`.
#' @return A list with a `vertices` field, containing the rownames of `X`, and then a list `graphs` one (named) entry for each value in `eps`. Each entry is a list with a `graph` field, storing the (undirected) edges in the Rips-Vietoris complex in matrix format, and a `clusters` field, containing vectors of the data indices (or row names) in each connected component of the Rips graph.
#' @export
#' @importFrom methods is
#' @importFrom stats dist
#' @importFrom foreach foreach %do%
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @examples
#'
#' if(require("TDA") & require("igraph"))
#' {
#'   # simulate data from the unit circle and calculate 
#'   # its diagram
#'   df <- TDA::circleUnif(n = 25)
#'   diag <- TDA::ripsDiag(df,maxdimension = 1,maxscale = 2)
#'   
#'   # get minimum death radius of any data cluster
#'   min_death_H0 <- 
#'   min(diag$diagram[which(diag$diagram[,1] == 0),3L])
#'   
#'   # get birth and death radius of the loop
#'   loop_birth <- as.numeric(diag$diagram[nrow(diag$diagram),2L])
#'   loop_death <- as.numeric(diag$diagram[nrow(diag$diagram),3L])
#'
#'   # compute Rips-Vietoris graphs at radii half of 
#'   # min_death_H0 and the mean of loop_birth and 
#'   # loop_death, returning clusters
#'   graphs <- rips_graphs(X = df,eps = 
#'   c(0.5*min_death_H0,(loop_birth + loop_death)/2))
#'
#'   # verify that there are 25 clusters for the smaller radius
#'   length(graphs$graphs[[0.5*min_death_H0]]$clusters)
#'   
#' }

rips_graphs <- function(X,distance_mat = FALSE,eps,return_clusters = TRUE){
  
  # function to compute static Rips complexes from a filtration at
  # specific epsilon radius values
  
  # error check parameters
  check_param(param = eps,param_name = "eps",numeric = T,whole_numbers = F,multiple = T,finite = T,positive = T)
  
  if(is.null(distance_mat))
  {
    stop("distance_mat must not be NULL.")
  }
  if(length(distance_mat) > 1 | !methods::is(distance_mat,"logical"))
  {
    stop("distance_mat must be a single logical (i.e. T or F).")
  }
  if(is.na(distance_mat) | is.nan(distance_mat) )
  {
    stop("distance_mat must not be NA/NAN.")
  }
  
  if(!methods::is(X,"data.frame") & !methods::is(X,"matrix"))
  {
    stop("X must either be a dataframe or a matrix.")
  }
  if(nrow(X) < 2 | ncol(X) < 1)
  {
    stop("X must have at least two rows and one column.")
  }
  if(length(which(stats::complete.cases(X) == F)) > 0)
  {
    stop("X must not contain any missing values.")
  }
  if(distance_mat == T & (ncol(X) != nrow(X) | !methods::is(X,"matrix")))
  {
    stop("if distance_mat is TRUE then X must be a square matrix.")
  }
  if((methods::is(X,"matrix") & !methods::is(X[1,1],"numeric")) | (methods::is(X,"data.frame") & length(which(unlist(lapply(X,is.numeric)))) < ncol(X)))
  {
    stop("X must have only numeric entries.")
  }
  
  if(is.null(return_clusters))
  {
    stop("return_clusters must not be NULL.")
  }
  if(length(return_clusters) > 1 | !methods::is(return_clusters,"logical"))
  {
    stop("return_clusters must be a single logical (i.e. T or F).")
  }
  if(is.na(return_clusters) | is.nan(return_clusters) )
  {
    stop("return_clusters must not be NA/NAN.")
  }
  
  # get rownames of X
  if(is.null(rownames(X)))
  {
    rowNames <- F
    vertices <- as.character(1:nrow(X))
  }else
  {
    rowNames <- T
    vertices <- rownames(X)
  }
  
  # if X is not a distance matrix, convert it to one
  if(!distance_mat)
  {
    X <- as.matrix(dist(X))
  }
  
  # for each value in eps, calculate complex (and optionally clusters)
  # in sequence
  ret_list <- foreach::`%do%`(obj = foreach::foreach(e = eps),ex = {
    
    # compute the complex defined by one particular eps value
    complex <- which(X <= e,arr.ind = T)
    
    # format and remove duplicate rows
    rownames(complex) <- NULL
    if(methods::is(complex,"integer"))
    {
      complex <- t(as.matrix(complex))
    }
    complex <- complex[which(complex[,1] < complex[,2]),]
    if(methods::is(complex,"integer"))
    {
      complex <- t(as.matrix(complex))
    }
    
    # replace row numbers with row names if possible
    if(rowNames & nrow(complex) > 0)
    {
      complex[,1] <- vertices[complex[,1]]
      complex[,2] <- vertices[as.numeric(complex[,2])]
    }
    
    ret_list <- list(graph = complex)
    
    # compute clusters if necessary
    if(return_clusters)
    {
      # initialize empty list of clusters
      clusters <- list()
      
      # perform DFS to find all clusters
      if(rowNames)
      {
        to_visit <- vertices
      }else
      {
        num_points <- nrow(X)
        to_visit <- 1:num_points
      }
      while(length(to_visit) > 0)
      {
        new_cluster <- c()
        to_visit_cluster <- c(to_visit[[1]])
        while(length(to_visit_cluster) > 0)
        {
          p <- to_visit_cluster[[1]]
          if(length(to_visit_cluster) > 1)
          {
            to_visit_cluster <- to_visit_cluster[2:length(to_visit_cluster)]
          }else
          {
            to_visit_cluster <- c()
          }
          new_cluster <- c(new_cluster,p)
          adj <- c(complex[which(complex[,1] == p),2L],complex[which(complex[,2] == p),1L])
          adj <- setdiff(adj,new_cluster)
          to_visit_cluster <- unique(c(to_visit_cluster,adj))
        }
        clusters[[length(clusters) + 1]] <- new_cluster
        to_visit <- setdiff(to_visit,new_cluster)
      }
      
      ret_list$clusters <- clusters
      
    }
    
    return(ret_list)
    
  })
  names(ret_list) <- as.character(eps)

  return(list(vertices = vertices,graphs = ret_list))
  
}

#### PLOT RIPS GRAPHS ####
#' Plot a Rips graph using the \link{igraph} package.
#' 
#' This function will throw an error if the \link{igraph} package is not installed.
#'
#' @param graphs the output of a `rips_graphs` function call.
#' @param eps the numeric radius of the graph in `graphs` to plot.
#' @param cols an optional character vector of vertex colors, default `NULL`.
#' @param component_of a vertex name, only the component of the graph containing that vertex will be plotted (useful for identifying representative (co)cycles in graphs). Default `NULL` (plot the whole graph).
#' @param plot_isolated_vertices a boolean representing whether or not to plot isolated vertices, default `FALSE`.
#' @param vertex_labels a boolean representing whether or not to plot vertex labels, default `TRUE`.
#' @param return_layout a boolean representing whether or not to return the plotting layout (x-y coordinates of each vertex) and the vertex labels, default `FALSE`.
#' @return if `return_layout` is `TRUE` then a list with elements "layout" (the numeric matrix of vertex x-y coordinates) and "vertices" (vertex labels), otherwise the function does not return anything.
#' @export
#' @importFrom methods is
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @examples
#'
#' if(require("TDA") & require("igraph"))
#' {
#'   # simulate data from the unit circle and calculate 
#'   # its diagram
#'   df <- TDA::circleUnif(n = 25)
#'   diag <- TDA::ripsDiag(df,maxdimension = 1,maxscale = 2)
#'   
#'   # get minimum death radius of any data cluster
#'   min_death_H0 <- 
#'   min(diag$diagram[which(diag$diagram[,1] == 0),3L])
#'   
#'   # get birth and death radius of the loop
#'   loop_birth <- as.numeric(diag$diagram[nrow(diag$diagram),2L])
#'   loop_death <- as.numeric(diag$diagram[nrow(diag$diagram),3L])
#'
#'   # compute Rips-Vietoris graphs at radii half of 
#'   # min_death_H0 and the mean of loop_birth and 
#'   # loop_death, returning clusters
#'   graphs <- rips_graphs(X = df,eps = 
#'   c(0.5*min_death_H0,(loop_birth + loop_death)/2))
#'   
#'   # plot graph of smaller (first) radius
#'   plot_rips_graph(graphs = graphs,eps = 0.5*min_death_H0)
#'   
#'   # plot graph of larger (second) radius
#'   plot_rips_graph(graphs = graphs,eps = (loop_birth + loop_death)/2)
#'   
#'   # repeat but with rownames for df, each vertex
#'   # will be plotted with its rownames
#'   rownames(df) <- paste0("V",1:25)
#'   graphs <- rips_graphs(X = df,eps = 
#'   c(0.5*min_death_H0,(loop_birth + loop_death)/2))
#'   plot_rips_graph(graphs = graphs,eps = 0.5*min_death_H0)
#'   
#'   # plot without vertex labels
#'   plot_rips_graph(graphs = graphs,eps = 0.5*min_death_H0,
#'                   vertex_labels = FALSE)
#'   
#'   # remove isolated vertices
#'   plot_rips_graph(graphs = graphs,eps = 0.5*min_death_H0,
#'                   plot_isolated_vertices = FALSE)
#'                  
#'   # plot only the graph component containing vertex "1"
#'   plot_rips_graph(graphs = graphs,eps = 0.5*min_death_H0,
#'                   component_of = 1)
#'  
#'   # return the layout of the graph for further
#'   # plotting customization
#'   plot_info <- plot_rips_graph(graphs = graphs,eps = 0.5*min_death_H0,
#'                                return_layout = TRUE)
#'   
#' }

plot_rips_graph <- function(graphs,eps,cols = NULL,component_of = NULL,plot_isolated_vertices = FALSE,return_layout = FALSE,vertex_labels = TRUE){
  
  # check for igraph
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package \"igraph\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  # error check parameters
  if(!is.list(graphs) | length(graphs) != 2 | length(which(names(graphs) == c("vertices","graphs"))) != 2 | !methods::is(graphs$graphs,"list"))
  {
    stop("graphs must be the output of a rips_graphs function call.")
  }
  check_param(param_name = "eps",param = eps,numeric = T,multiple = F,finite = T,positive = T)
  if(eps %in% names(graphs$graphs) == F)
  {
    stop("eps must be the scale of one of the graphs.")
  }
  
  if(!is.null(cols))
  {
    check_param(param_name = "cols",param = cols,multiple = T)
    if(!is.character(cols))
    {
      stop("cols must be a character vector.")
    }
    if(length(which(is.na(cols) | is.nan(cols) | is.null(cols))) > 0)
    {
      stop("cols must not having missing values.")
    }
    if(length(cols) != length(graphs$vertices))
    {
      stop("cols must have the same length as graphs$vertices.")
    }
  }
  
  if(is.null(plot_isolated_vertices))
  {
    stop("plot_isolated_vertices must not be NULL.")
  }
  if(length(plot_isolated_vertices) > 1 | !methods::is(plot_isolated_vertices,"logical"))
  {
    stop("plot_isolated_vertices must be a single boolean value.")
  }
  if(is.na(plot_isolated_vertices) | is.nan(plot_isolated_vertices) )
  {
    stop("plot_isolated_vertices must not be NA/NAN.")
  }
  
  if(is.null(return_layout))
  {
    stop("return_layout must not be NULL.")
  }
  if(length(return_layout) > 1 | !methods::is(return_layout,"logical"))
  {
    stop("return_layout must be a single boolean value.")
  }
  if(is.na(return_layout) | is.nan(return_layout) )
  {
    stop("return_layout must not be NA/NAN.")
  }
  
  if(is.null(vertex_labels))
  {
    stop("vertex_labels must not be NULL.")
  }
  if(length(vertex_labels) > 1 | !methods::is(vertex_labels,"logical"))
  {
    stop("vertex_labels must be a single boolean value.")
  }
  if(is.na(vertex_labels) | is.nan(vertex_labels) )
  {
    stop("vertex_labels must not be NA/NAN.")
  }
  
  if(!is.null(component_of))
  {
    check_param(param_name = "component_of",param = component_of,multiple = F)
    if(is.na(component_of) | is.nan(component_of))
    {
      stop("component_of must not be NA or NaN.")
    }
    if(component_of %in% graphs$vertices == F)
    {
      stop("component_of must be an element of graphs$vertices.")
    }
  }
  
  # create graph
  g <- igraph::graph_from_data_frame(graphs$graphs[[which(names(graphs$graphs) == eps)]]$graph,directed = FALSE,vertices = graphs$vertices)
  if(!is.null(cols))
  {
    igraph::V(g)$color <- cols
  }
  
  # if desired, subset by component of one vertex
  if(!is.null(component_of))
  {
    components <- igraph::components(g)
    comp <- graphs$vertices[which(components$membership == components$membership[[component_of]])]
    g <- igraph::induced_subgraph(g,comp)
  }
  
  # if desired remove isolated vertices
  if(plot_isolated_vertices == F)
  {
    isolated = which(igraph::degree(g)==0)
    g = igraph::delete.vertices(g, isolated)
    if(length(igraph::V(g)) == 0)
    {
      warning("graph had 0 non-isolated vertices, empty graph will be plotted.")
    }
  }
  
  # get graph layout
  layout <- igraph::layout_nicely(g)
  layout <- apply(X = layout,MARGIN = 2L,FUN = function(X){
    
    return(-1 + 2*(X - min(X))/(max(X) - min(X)))
    
  })
  
  # plot graph
  if(vertex_labels)
  {
    igraph::plot.igraph(g,layout = layout)
  }else
  {
    igraph::plot.igraph(g,layout = layout,vertex.label = NA)
  }
  
  # if desired return layout and final vertex labels
  if(return_layout)
  {
    return(list(layout = layout,vertices = igraph::V(g)))
  }
  
}
