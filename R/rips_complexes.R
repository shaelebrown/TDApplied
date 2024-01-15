#### COMPUTE STATIC SIMPLICIAL COMPLEXES WITHIN RIPS FILTRATION ####
#' Compute Vietoris-Rips graphs of a dataset at particular epsilon radius values.
#'
#' Persistence diagrams computed from Rips-Vietoris filtrations contain information about 
#' distance radius scales at which topological features of a dataset exist, but the features
#' can be challenging to visualize, analyze and interpret. In order to help solve this problem the `vr_graphs`
#' function computes the 1-skeleton (i.e. graph) of Rips complexes at particular radii, called "Vietoris-Rips graphs" (VR graphs) in the literature.
#' 
#' This function may be used in conjunction with the igraph package to visualize the graphs (see \code{\link{plot_vr_graph}}).
#'
#' @param X either a point cloud data frame/matrix, or a distance matrix.
#' @param distance_mat a boolean representing if the input `X` is a distance matrix, default value is `FALSE`.
#' @param eps a numeric vector of the positive scales at which to compute the Rips-Vietoris complexes, i.e. all edges at most the specified values.
#' @param return_clusters a boolean determining if the connected components (i.e. data clusters) of the complex should be explicitly returned, default is `TRUE`.
#' @return A list with a `vertices` field, containing the rownames of `X`, and then a list `graphs` one (named) entry for each value in `eps`. Each entry is a list with a `graph` field, storing the (undirected) edges in the Rips-Vietoris complex in matrix format, and a `clusters` field, containing vectors of the data indices (or row names) in each connected component of the Rips graph.
#' @export
#' @importFrom stats dist
#' @importFrom foreach foreach %do%
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @references 
#' 
#' A Zomorodian, The tidy set: A minimal simplicial set for computing homology of clique complexes in Proceedings of the Twenty-Sixth Annual Symposium on Computational Geometry, SoCG ’10. (Association for Computing Machinery, New York, NY, USA), p. 257–266 (2010).
#' 
#' @examples
#'
#' if(require("TDAstats") & require("igraph"))
#' {
#'   # simulate data from the unit circle and calculate 
#'   # its diagram
#'   df <- TDAstats::circle2d[sample(1:100,25),]
#'   diag <- TDAstats::calculate_homology(df,
#'                                        dim = 1,
#'                                        threshold = 2)
#'   
#'   # get minimum death radius of any data cluster
#'   min_death_H0 <- 
#'   min(diag[which(diag[,1] == 0),3L])
#'   
#'   # get birth and death radius of the loop
#'   loop_birth <- as.numeric(diag[nrow(diag),2L])
#'   loop_death <- as.numeric(diag[nrow(diag),3L])
#'
#'   # compute VR graphs at radii half of 
#'   # min_death_H0 and the mean of loop_birth and 
#'   # loop_death, returning clusters
#'   graphs <- vr_graphs(X = df,eps = 
#'   c(0.5*min_death_H0,(loop_birth + loop_death)/2))
#'
#'   # verify that there are 25 clusters for the smaller radius
#'   length(graphs$graphs[[1]]$clusters)
#'   
#' }

vr_graphs <- function(X,distance_mat = FALSE,eps,return_clusters = TRUE){
  
  # function to compute static Rips complexes from a filtration at
  # specific epsilon radius values
  
  # to avoid build issues
  e <- NULL
  
  # error check parameters
  check_param(param = eps,param_name = "eps",numeric = T,whole_numbers = F,multiple = T,finite = T,positive = T)
  
  if(is.null(distance_mat))
  {
    stop("distance_mat must not be NULL.")
  }
  if(length(distance_mat) > 1 | !inherits(distance_mat,"logical"))
  {
    stop("distance_mat must be a single logical (i.e. T or F).")
  }
  if(is.na(distance_mat) | is.nan(distance_mat) )
  {
    stop("distance_mat must not be NA/NAN.")
  }
  
  if(!inherits(X,"data.frame") & !inherits(X,"matrix"))
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
  if(distance_mat == T & (ncol(X) != nrow(X) | !inherits(X,"matrix")))
  {
    stop("if distance_mat is TRUE then X must be a square matrix.")
  }
  if((inherits(X,"matrix") & !inherits(X[1,1],"numeric")) | (inherits(X,"data.frame") & length(which(unlist(lapply(X,is.numeric)))) < ncol(X)))
  {
    stop("X must have only numeric entries.")
  }
  
  if(is.null(return_clusters))
  {
    stop("return_clusters must not be NULL.")
  }
  if(length(return_clusters) > 1 | !inherits(return_clusters,"logical"))
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
    if(inherits(complex,"integer"))
    {
      complex <- t(as.matrix(complex))
    }
    complex <- complex[which(complex[,1] < complex[,2]),]
    if(inherits(complex,"integer"))
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
#' Plot a VR graph using the igraph package.
#' 
#' This function will throw an error if the igraph package is not installed.
#'
#' @param graphs the output of a `vr_graphs` function call.
#' @param eps the numeric radius of the graph in `graphs` to plot.
#' @param cols an optional character vector of vertex colors, default `NULL`.
#' @param layout an optional 2D matrix of vertex coordinates, default `NULL`. If row names are supplied they can be used to subset a graph by those vertex names.
#' @param title an optional str title for the plot, default `NULL`.
#' @param component_of a vertex name (integer or character), only the component of the graph containing that vertex will be plotted (useful for identifying representative (co)cycles in graphs). Default `NULL` (plot the whole graph).
#' @param plot_isolated_vertices a boolean representing whether or not to plot isolated vertices, default `FALSE`.
#' @param vertex_labels a boolean representing whether or not to plot vertex labels, default `TRUE`.
#' @param return_layout a boolean representing whether or not to return the plotting layout (x-y coordinates of each vertex) and the vertex labels, default `FALSE`.
#' @return if `return_layout` is `TRUE` then a list with elements "layout" (the numeric matrix of vertex x-y coordinates) and "vertices" (character vertex labels), otherwise the function does not return anything.
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @examples
#'
#' if(require("TDAstats") & require("igraph"))
#' {
#'   # simulate data from the unit circle and calculate 
#'   # its diagram
#'   df <- TDAstats::circle2d[sample(1:100,25),]
#'   diag <- TDAstats::calculate_homology(df,
#'                                        dim = 1,
#'                                        threshold = 2)
#'   
#'   # get minimum death radius of any data cluster
#'   min_death_H0 <- 
#'   min(diag[which(diag[,1] == 0),3L])
#'   
#'   # get birth and death radius of the loop
#'   loop_birth <- as.numeric(diag[nrow(diag),2L])
#'   loop_death <- as.numeric(diag[nrow(diag),3L])
#'
#'   # compute VR graphs at radii half of 
#'   # min_death_H0 and the mean of loop_birth and 
#'   # loop_death, returning clusters
#'   graphs <- vr_graphs(X = df,eps = 
#'   c(0.5*min_death_H0,(loop_birth + loop_death)/2))
#'   
#'   # plot graph of smaller (first) radius
#'   plot_vr_graph(graphs = graphs,eps = 0.5*min_death_H0,
#'                   plot_isolated_vertices = TRUE)
#'   
#'   # plot graph of larger (second) radius
#'   plot_vr_graph(graphs = graphs,eps = (loop_birth + loop_death)/2)
#'   
#'   # repeat but with rownames for df, each vertex
#'   # will be plotted with its rownames
#'   rownames(df) <- paste0("V",1:25)
#'   graphs <- vr_graphs(X = df,eps = 
#'   c(0.5*min_death_H0,(loop_birth + loop_death)/2))
#'   plot_vr_graph(graphs = graphs,eps = 0.5*min_death_H0,
#'                   plot_isolated_vertices = TRUE)
#'   
#'   # plot without vertex labels
#'   plot_vr_graph(graphs = graphs,eps = (loop_birth + loop_death)/2,
#'                   vertex_labels = FALSE)
#'                  
#'   # plot only the graph component containing vertex "1"
#'   plot_vr_graph(graphs = graphs,eps = 0.5*min_death_H0,
#'                   component_of = "V1",plot_isolated_vertices = TRUE)
#'  
#'   # save the layout of the graph for adding features to 
#'   # the same graph layout, like color
#'   layout <- plot_vr_graph(graphs = graphs,eps = (loop_birth + loop_death)/2,
#'                             return_layout = TRUE,vertex_labels = TRUE)
#'   cols <- rep("blue",25)
#'   cols[1:5] <- "red"
#'   plot_vr_graph(graphs = graphs,eps = (loop_birth + loop_death)/2,cols = cols,
#'                   layout = layout)
#'   
#' }

plot_vr_graph <- function(graphs,eps,cols = NULL,layout = NULL,title = NULL,component_of = NULL,plot_isolated_vertices = FALSE,return_layout = FALSE,vertex_labels = TRUE){
  
  # check for igraph
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "Package \"igraph\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
  # error check parameters
  if(!is.list(graphs) | length(graphs) != 2 | length(which(names(graphs) == c("vertices","graphs"))) != 2 | !inherits(graphs$graphs,"list"))
  {
    stop("graphs must be the output of a vr_graphs function call.")
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
  
  if(!is.null(layout))
  {
    if(!is.matrix(layout))
    {
      stop("layout must be a matrix.")
    }
    if(ncol(layout) != 2)
    {
      stop("layout must have 2 columns.")
    }
    if(is.null(rownames(layout)))
    {
      if(nrow(layout) != length(graphs$vertices))
      {
        stop("when layout has no rownames, layout must have the same number of rows as the number of vertices in the graph.")
      }
    }else
    {
      # subset graphs by row names of layout
      rnames <- as.character(rownames(layout))
      inds <- unlist(lapply(X = rnames,FUN = function(X){
        
        ind <- which(graphs$vertices == X)
        if(length(ind) == 0)
        {
          return(NA)
        }
        return(ind)
        
      }))
      layout <- layout[which(!is.na(inds)),]
      if(is.null(dim(layout)))
      {
        layout <- matrix(data = layout,ncol = 2)
        rownames(layout) <- rnames[which(!is.na(inds))]
      }
      cols <- cols[which(!is.na(inds))]
      inds <- inds[which(!is.na(inds))]
      graphs$graphs <- lapply(X = graphs$graphs,FUN = function(X){
        
        return(list(graph = X$graph[which(as.character(X$graph[,1]) %in% rnames & as.character(X$graph[,2]) %in% rnames),]))
        
      })
      graphs$vertices <- rnames
    }
  }
  
  if(!is.null(title))
  {
    if(!is.character(title) | length(title) > 1)
    {
      stop("title must be a single character string.")
    }
  }
  
  if(is.null(plot_isolated_vertices))
  {
    stop("plot_isolated_vertices must not be NULL.")
  }
  if(length(plot_isolated_vertices) > 1 | !inherits(plot_isolated_vertices,"logical"))
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
  if(length(return_layout) > 1 | !inherits(return_layout,"logical"))
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
  if(length(vertex_labels) > 1 | !inherits(vertex_labels,"logical"))
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
    if(as.character(component_of) %in% graphs$vertices == F)
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
      warning("graph had 0 non-isolated vertices, empty graph will be plotted and no layout will be returned.")
      if(is.null(title))
      {
        if(vertex_labels)
        {
          igraph::plot.igraph(g,margin = 0)
        }else
        {
          igraph::plot.igraph(g,vertex.label = NA,margin = 0)
        }
      }else
      {
        if(vertex_labels)
        {
          igraph::plot.igraph(g,margin = 0,main = title)
        }else
        {
          igraph::plot.igraph(g,vertex.label = NA,margin = 0,main = title)
        }
      }
    }
  }
  
  # get graph layout if needed and standardize
  if(is.null(layout))
  {
    layout <- igraph::layout_nicely(g)
  }else
  {
    layout <- layout[which(rownames(layout) %in% names(igraph::V(g))),]
    if(is.null(dim(layout)))
    {
      layout <- matrix(data = layout,nrow = 1,ncol = 2)
    }
  }
  if(nrow(layout) > 1)
  {
    layout <- apply(X = layout,MARGIN = 2L,FUN = function(X){
      
      return(-1 + 2*(X - min(X))/(max(X) - min(X)))
      
    }) 
  }
  
  # plot graph
  if(is.null(title))
  {
    if(vertex_labels)
    {
      igraph::plot.igraph(g,layout = layout,margin = 0)
    }else
    {
      igraph::plot.igraph(g,layout = layout,vertex.label = NA,margin = 0)
    }
  }else
  {
    if(vertex_labels)
    {
      igraph::plot.igraph(g,layout = layout,margin = 0,main = title)
    }else
    {
      igraph::plot.igraph(g,layout = layout,vertex.label = NA,margin = 0,main = title)
    }
  }
  
  
  # if desired return layout and final vertex labels
  if(return_layout)
  {
    rownames(layout) <- names(igraph::V(g))
    return(layout)
  }
  
}
