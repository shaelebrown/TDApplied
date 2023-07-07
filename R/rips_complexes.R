#### COMPUTE STATIC SIMPLICIAL COMPLEXES WITHIN RIPS FILTRATION ####
#' Compute Rips graphs of a dataset at particular epsilon radius values.
#'
#' Persistence diagrams computed from Rips-Vietoris filtrations contain information about 
#' distance radius scales at which topological features of a dataset exist, but the features
#' can be challenging to visualize, analyze and interpret. In order to help solve this problem the `rips_graphs`
#' function computes the 1-skeleton (i.e. graph) of Rips complexes at particular radii.
#' 
#' This function may be used in conjunction with the \link{igraph} package to visualize the graphs (see examples).
#'
#' @param X either a point cloud data frame/matrix, or a distance matrix.
#' @param distance_mat a boolean representing if the input `X` is a distance matrix, default value is `FALSE`.
#' @param eps a numeric vector of the positive scales at which to compute the Rips-Vietoris complexes, i.e. all edges at most the specified values.
#' @param return_clusters a boolean determining if the connected components (i.e. data clusters) of the complex should be explicitly returned, default is `TRUE`.
#' @return A list with a `data` field, containing a dataframe with the rownames of `X`, and then a list `graphs` one entry for each value in `eps`. Each entry is a list with a `graph` field, storing the (undirected) edges in the Rips-Vietoris complex in matrix format, and a `clusters` field, containing vectors of the data indices in each connected component of the Rips-Vietoris complex.
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
#'   graph <- rips_graphs(X = df,eps = 
#'   c(0.5*min_death_H0,(loop_birth + loop_death)/2))
#'
#'   # verify that there are 25 clusters for the smaller radius
#'   length(graph$graphs[[1]]$clusters)
#'   
#'   # plot graph of larger radius using 
#'   # the igraph package to see the loop
#'   g <- igraph::graph_from_data_frame(
#'   graph$graphs[[2]]$graph,
#'   directed = FALSE,vertices = graph$data)
#'   igraph::plot.igraph(g)
#'   
#' }

rips_graphs <- function(X,distance_mat = FALSE,eps,return_clusters = TRUE){
  
  # function to compute static Rips complexes from a filtration at
  # specific epsilon radius values
  
  # error check parameters
  check_param(param = eps,param_name = "eps",numeric = T,whole_numbers = F,multiple = T,finite = T,non_negative = T,positive = T)
  
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
    data <- data.frame(name = as.character(1:nrow(X)))
  }else
  {
    data <- data.frame(name = rownames(X))
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
    
    ret_list <- list(graph = complex)
    
    # compute clusters if necessary
    if(return_clusters)
    {
      # initialize empty list of clusters
      clusters <- list()
      
      # perform DFS to find all clusters
      num_points <- nrow(X)
      to_visit <- 1:num_points
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

  return(list(data = data,graphs = ret_list))
  
}
