#### DIAGRAM DISTANCE METRICS ####
#' Calculate distances between pairs of persistence diagrams
#'
#' Calculates the distance between a pair of persistce diagrams in a particular homological dimension,
#' either the output from a \code{\link{diagram_to_df}} function call or from a TDA homology calculation like ripsDiag.
#' Different TDA sources define distances
#' differently, and this function has functionality to compute distances like
#' in the R package TDA (based on the C++ library Dionysus, see
#' \url{https://mrzv.org/software/dionysus2/}) or like in the
#' paper for kernel calculations of persistence diagrams (
#' \url{https://proceedings.neurips.cc/paper/2018/file/959ab9a0695c467e7caf75431a872e5c-Paper.pdf}).
#'
#' The `D1` and `D2` parameters should be persistence diagrams, outputted
#' from a homology calculation in the package TDA, or such a
#' persistence diagram converted to a data frame via the function \code{\link{diagram_to_df}}.
#' The `dim` parameter should be a positive finite integer.
#' The `p` parameter should be a positive integer or Inf. The `distance` parameter
#' should be a string, either "wasserstein" or "fisher". The `sigma` parameter is a single positive number representing the bandwith for the Fisher information metric.
#'
#' @param D1 the first persistence diagram, either computed from TDA or converted to a data frame with diagram_to_df.
#' @param D2 the second persistence diagram, either computed from TDA or converted to a data frame with diagram_to_df.
#' @param dim the homological dimension in which the distance is to be computed.
#' @param p  the wasserstein power parameter. Default value is 2.
#' @param distance a string which determines which type of distance calculation to carry out, either "wasserstein" (default) or "fisher".
#' @param sigma either NULL (default) or a positive number representing the bandwith for the Fisher information metric
#'
#' @return the numeric value of the distance calculation.
#' @importFrom rdist cdist
#' @importFrom clue solve_LSAP
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @references
#' Robinson T, Turner K (2017). "Hypothesis testing for topological data analysis." \url{https://link.springer.com/article/10.1007/s41468-017-0008-7}.
#' 
#' Le T, Yamada M (2018). "Persistence fisher kernel: a riemannian manifold kernel for persistence diagrams." \url{https://proceedings.neurips.cc/paper/2018/file/959ab9a0695c467e7caf75431a872e5c-Paper.pdf}.
#' @examples
#'
#' # create two diagrams with package TDA based on 2D Gaussians
#' diag1 <- TDA::ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),
#' y = rnorm(100,mean = 0,sd = 1)),
#' maxscale = 1,
#' maxdimension = 1)
#' diag2 <- TDA::ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),
#' y = rnorm(100,mean = 0,sd = 1)),
#' maxscale = 1,
#' maxdimension = 1)
#'
#' # calculate their wasserstein distance
#' wass <- diagram_distance(D1 = diag1,D2 = diag2,dim = 1,p = 2,distance = "wasserstein")
#'
#' # calculate their bottleneck distance
#' bottleneck <- diagram_distance(D1 = diag2,D2 = diag2,dim = 1,p = Inf,distance = "wasserstein")
#'
#' # repeat wasserstein calculation but with diagrams converted to data frames
#' diag1_df <- diagram_to_df(d = diag1)
#' diag2_df <- diagram_to_df(d = diag2)
#' wass_df <- diagram_distance(D1 = diag1_df,D2 = diag2_df,dim = 1,p = 2,distance = "wasserstein")
#' 
#' # now do Fisher information metric calculation
#' fisher_df <- diagram_distance(D1 = diag1_df,D2 = diag2_df,dim = 1,distance = "fisher",sigma = 1)

diagram_distance <- function(D1,D2,dim,p = 2,distance = "wasserstein",sigma = 1){

  # function to compute the wasserstein/bottleneck/Fisher information metric between two diagrams
  # D1 and D2 are diagrams, possibly stored as data frames
  # dim is the dimension to subset
  # p is the power of the wasserstein distance, p >= 1
  # distance is either "wasserstein" (default) or "fisher"
  # sigma is the positive bandwidth for the Fisher information metric, default 1

  # for standalone usage force D1 and D2 to be data frames if they are the output of a homology calculation
  if(is.list(D1) && length(D1) == 1 && names(D1) == "diagram" && class(D1$diagram) == "diagram")
  {
    # D1 is the output from a TDA calculation
    D1 <- diagram_to_df(D1)
  }else
  {
    if(class(D1) != "data.frame")
    {
      stop("D1 must either be the output of a TDA computation or data frame.")
    }
  }

  # error check D1
  check_diagram(D1)

  if(is.list(D2) && length(D2) == 1 && names(D2) == "diagram" && class(D2$diagram) == "diagram")
  {
    # D1 is the output from a TDA calculation
    D2 <- diagram_to_df(D2)
  }else
  {
    if(class(D2) != "data.frame")
    {
      stop("D2 must either be the output of a TDA computation or data frame.")
    }
  }

  # error check D2
  check_diagram(D2)
  
  # error check dim
  if(is.null(dim) | is.na(dim) | is.nan(dim) | !is.numeric(dim) | length(dim) > 1 | as.integer(dim) != dim | dim < 0)
  {
    stop("dim must be a non-negative whole number.")
  }
  
  # error check distance
  if(is.null(distance) | !is.character(distance) | length(distance) > 1 | distance %in% c("wasserstein","fisher") == F)
  {
    stop("distance must be either \'wasserstein\' or \'fisher\'.")
  }
  
  # error check p
  if((is.null(p) & distance == "wasserstein") | is.na(p) | is.nan(p) | !is.numeric(p) | length(p) > 1 | p < 1)
  {
    stop("p must be a number at least one for the wasserstein metric.")
  }
  
  # if persistence Fisher distance, sigma must be a positive number
  if(distance == "fisher" & !is.null(sigma))
  {
    if(is.na(sigma) | is.nan(sigma) | !is.numeric(sigma) | length(sigma) > 1 | sigma <= 0)
    {
      stop("For persistence Fisher distance sigma must be a positive number.") 
    }
  }

  # subset both diagrams by dimension dim and for birth and death columns
  D1_subset <- D1[which(D1$dimension == dim),1:3]
  D2_subset <- D2[which(D2$dimension == dim),1:3]
  D1_subset <- D1_subset[,2:3]
  D2_subset <- D2_subset[,2:3]

  # create empty diagonals for the persistence diagrams
  diag1 <- D1_subset[0,]
  diag2 <- D2_subset[0,]

  # if both subsets are empty then set their distance to 0
  if(nrow(D1_subset) == 0 & nrow(D2_subset) == 0)
  {
    return(0)
  }

  # remove diagonal entries from D1_subset and D2_subset
  D1_subset <- D1_subset[which(D1_subset[,1] != D1_subset[,2]),]
  D2_subset <- D2_subset[which(D2_subset[,1] != D2_subset[,2]),]

  # for each non-trivial element in D1_subset we add its projection onto the diagonal in diag1
  if(nrow(D1_subset) > 0)
  {
    for(i in 1:nrow(D1_subset))
    {
      proj_diag <- mean(as.numeric(D1_subset[i,]))
      diag1 <- rbind(diag1,data.frame(birth = proj_diag,death = proj_diag))
    }
  }

  # for each non-trivial element in D2_subset we add its projection onto the diagonal in diag2
  if(nrow(D2_subset) > 0)
  {
    for(i in 1:nrow(D2_subset))
    {
      proj_diag <- mean(as.numeric(D2_subset[i,]))
      diag2 <- rbind(diag2,data.frame(birth = proj_diag,death = proj_diag))
    }
  }

  # since an element b of D1_subset is either matched to an element of D2 or to the projection of b onto the diagonal
  # we form the two sets to be matched by row binding D1_subset with diag2 and D2_subset with diag1
  D1_subset <- rbind(D1_subset,diag2)
  D2_subset <- rbind(D2_subset,diag1)

  if(distance == "wasserstein")
  {
    # compute the bottleneck distance matrix between rows of D1_subset and D2_subset
    dist_mat <- as.matrix(rdist::cdist(D1_subset,D2_subset,metric = "maximum"))
    
    # use the Hungarian algorithm from the clue package to find the minimal weight matching
    best_match <- as.numeric(clue::solve_LSAP(x = dist_mat,maximum = F))
    seq_match <- 1:length(best_match)
    
    # subset best match by removing all pairs between diagonal points
    indices <- cbind(seq_match,best_match)
    indices <- indices[which(indices[,1] <= (nrow(D1_subset) - nrow(diag2)) | indices[,2] <= (nrow(D2_subset) - nrow(diag1))),]
    
    # if p is finite, exponentiate each matched distance and take the p-th root of their sum
    if(is.finite(p))
    {
      return((sum(dist_mat[indices]^(p)))^(1/p))
    }
    
    # otherwise, return the regular bottleneck distance
    return(max(dist_mat[indices]))
    
  }else
  {
    # compute the persistence Fisher distance
    
    # get all unique points in both diagrams and their diagonal projections
    theta <- unique(rbind(D1_subset,D2_subset))
    
    rho_1 <- unlist(lapply(X = 1:nrow(theta),FUN = function(X){
      
      x <- as.numeric(theta[X,])
      sum <- 0
      for(i in 1:nrow(D1_subset))
      {
        # compute mvn normal pdf values and sum up
        u <- as.numeric(D1_subset[i,])
        sum <- sum + exp(-1*((x[[1]]-u[[1]])^2 + (x[[2]]-u[[2]])^2)/(2*sigma^2))/(2*pi*sigma^2)
      }
      return(sum)
      
    }))
    rho_1 <- rho_1/sum(rho_1)
    
    rho_2 <- unlist(lapply(X = 1:nrow(theta),FUN = function(X){
      
      x <- as.numeric(theta[X,])
      sum <- 0
      for(i in 1:nrow(D2_subset))
      {
        # compute mvn normal pdf values and sum up
        u <- as.numeric(D2_subset[i,])
        sum <- sum + exp(-1*((x[[1]]-u[[1]])^2 + (x[[2]]-u[[2]])^2)/(2*sigma^2))/(2*pi*sigma^2)
      }
      return(sum)
      
    }))
    rho_2 <- rho_2/sum(rho_2)
    
    # return dot product of elementwise square root
    return(as.numeric(sqrt(rho_1) %*% sqrt(rho_2)))
    
  }

}

#### DISTANCE MATRIX ####
#' Compute a distance matrix between persistence diagrams
#'
#' Calculate the distance matrix d for either a single list of persistence diagrams \eqn{D = (D_1,D_2,\dots,D_n)}, i.e. \eqn{d[i,j] = d(D_i,D_j)}, 
#' or between two lists, \eqn{D = (D_1,D_2,\dots,D_n)} and \eqn{D' = (D'_1,D'_2,\dots,D'_n)}, \eqn{d[i,j] = d(D_i,D'_j)}, in parallel.
#'
#' `diagrams` is a list of persistence diagrams and `other_diagrams` is the optional second list of persistence diagrams.
#' The `dim` parameter should be a positive finite integer, being the homological dimension in which to compute distances.
#' The `distance` parameter is the string determining which distance metric to use, `p` is the 
#' wasserstein power parameter, and
#' `t` is the positive scale parameter for the persistence Fisher kernel.
#'
#' @param diagrams the list of persistence diagrams, either the output from TDA calculations or the \code{\link{diagram_to_df}} function.
#' @param other_diagrams either NULL (default) or another list of persistence diagrams to compute a cross-distance matrix.
#' @param dim the homological dimension in which the distance is to be computed.
#' @param distance a character determining which metric to use, either "wasserstein" (default) or "fisher".
#' @param p the positive wasserstein power, default 2.
#' @param sigma a positive number representing the bandwith of the Fisher information metric, default NULL.
#'
#' @return the numeric distance matrix.
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ
#' @importFrom parallelly availableCores
#' @importFrom doParallel registerDoParallel
#' @importFrom iterators iter
#' @examples
#'
#' # create ten diagrams with package TDA based on 2D Gaussians
#' g <- lapply(X = 1:10,FUN = function(X){
#'
#' diag <- TDA::ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),
#' y = rnorm(100,mean = 0,sd = 1)),
#' maxscale = 1,
#' maxdimension = 1)
#' df <- diagram_to_df(d = diag)
#' return(df)
#'
#' })
#'
#' # calculate their distance matrix in dimension 1 with the 2-wasserstein metric
#' D <- distance_matrix(diagrams = g,dim = 1,distance = "wasserstein",p = 2)
#' 
#' # now do the cross distance matrix, should be the same as the original
#' D_cross <- distance_matrix(diagrams = g,other_diagrams = g,dim = 1,distance = "wasserstein",p = 2)

distance_matrix <- function(diagrams,other_diagrams = NULL,dim = 0,distance = "wasserstein",p = 2,sigma = NULL){
  
  # set internal variables to NULL to avoid build issues
  r <- NULL
  X <- NULL
  
  # error check diagrams argument
  if(is.null(diagrams))
  {
    stop("diagrams must be a list of persistence diagrams.")
  }
  if(!is.list(diagrams) | length(diagrams) < 2)
  {
    stop("diagrams must be a list of persistence diagrams of length at least 2.")
  }
  diagrams <- all_diagrams(diagram_groups = list(diagrams),inference = "independence")[[1]]
  
  # error check other_diagrams argument
  if(!is.null(other_diagrams))
  {
    if(!is.list(other_diagrams) | length(other_diagrams) < 2)
    {
      stop("diagrams must be a list of persistence diagrams of length at least 2.")
    }
    other_diagrams <- all_diagrams(diagram_groups = list(other_diagrams),inference = "independence")[[1]]
  }
  
  # check other parameters
  check_params(iterations = 10,p = p,q = 2,dims = c(dim),paired = F,distance = distance,sigma = sigma)
  
  # compute distance matrix in parallel
  m = length(diagrams)
  num_workers <- parallelly::availableCores(omit = 1)
  cl <- parallel::makeCluster(num_workers)
  doParallel::registerDoParallel(cl)
  parallel::clusterEvalQ(cl,c(library(clue),library(rdist)))
  parallel::clusterExport(cl,c("diagram_distance"))
  force(diagrams) # required for parallel computation in this environment
  
  if(is.null(other_diagrams))
  {
    d <- matrix(data = 0,nrow = m,ncol = m)
    d_off_diag <- foreach::`%dopar%`(obj = foreach::foreach(r = iterators::iter(which(upper.tri(d),arr.ind = T),by = 'row'),.combine = c),ex = {
      
      return(diagram_distance(D1 = diagrams[[r[[1]]]],D2 = diagrams[[r[[2]]]],dim = dim,p = p,distance = distance,sigma = sigma))
      
    })
    d[upper.tri(d)] <- d_off_diag
    d[which(upper.tri(d),arr.ind = T)[,c("col","row")]] <- d_off_diag
  }else
  {
    if(length(other_diagrams) > length(diagrams))
    {
      d <- foreach::`%dopar%`(obj = foreach::foreach(r = 1:length(other_diagrams),.combine = cbind),ex = {
        
        return(unlist(lapply(X = 1:length(diagrams),FUN = function(X){return(diagram_distance(D1 = other_diagrams[[r]],D2 = diagrams[[X]],dim = dim,p = p,distance = distance,sigma = sigma))})))
        
      })
    }else
    {
      d <- foreach::`%do%`(obj = foreach::foreach(r = 1:length(other_diagrams),.combine = cbind),ex = {
        
        return(foreach::`%dopar%`(obj = foreach::foreach(X = 1:length(diagrams),.combine = c),ex = {
          
          return(diagram_distance(D1 = other_diagrams[[r]],D2 = diagrams[[X]],dim = dim,p = p,distance = distance,sigma = sigma))
          
        }))
        
      })
    }
  }
  
  
  parallel::stopCluster(cl)
  
  return(d)
  
}

#### LOSS FUNCTION FOR GROUPS OF PERSISTENCE DIAGRAMS ####
#' Turner loss function for a set of groups of persistence diagrams
#'
#' Calculate the normalized sum of within-group exponentiated distances between pairs of persistence diagrams (stored as data frames)
#' for an arbitrary number of groups in parallel. The loss function is described in Robinson and Turner 2017
#' \url{https://link.springer.com/article/10.1007/s41468-017-0008-7}, but mathematically we
#' compute the distances of each within-group pair of persistence diagrams and exponentiate
#' each distance by \eqn{q} and take the \eqn{q}-th root of the sum.
#'
#' The `diagram_groups` parameter should be a list or vector of persistence diagrams stored as data frames.
#' The `dims` parameter is a vector of non-negative whole numbers, representing the homological dimensions
#' to calculate the loss function in. The `dist_mats` parameter should be a list, with one element for each element in the parameter `dims`,
#' which stores a matrix of distance calculations (with -1 entries for distance calculations yet to be completed).
#' The `p` parameter should be a number at least 1 and possibly Inf.
#' The `q` parameter should be a finite number at least 1. The `distance` parameter should be a string
#' either "wasserstein" or "fisher". The `sigma` parameter is the positive bandwith for the persistence
#' Fisher distance.
#'
#' @param diagram_groups groups (lists/vectors) of persistence diagrams, stored as lists of a data frame and
#'                          an index of the diagram in all the diagrams across all groups.
#' @param dist_mats distance matrices between all possible pairs of persistence diagrams across and within groups
#'                      storing the current distances which have been precomputed.
#' @param dims a numeric vector of which homological dimensions in which the loss function is to be computed.
#' @param p a positive wasserstein parameter, if Inf then the bottleneck distance.
#' @param q a finite exponent at least 1.
#' @param distance a string which determines which type of distance calculation to carry out, either "wasserstein" (default) or "fisher".
#' @param sigma the positive bandwith for the persistence Fisher distance.
#'
#' @importFrom parallel makeCluster clusterEvalQ clusterExport stopCluster
#' @importFrom parallelly availableCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom utils combn
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @references
#' Robinson T, Turner K (2017). "Hypothesis testing for topological data analysis." \url{https://link.springer.com/article/10.1007/s41468-017-0008-7}.
#' @return the numeric value of the Turner loss function.

loss <- function(diagram_groups,dist_mats,dims,p,q,distance,sigma){

  # function to compute the F_{p,q} loss between groups of diagrams
  # diagram_groups are the (possibly permuted) groups of diagrams
  # dist mats stores the current distance matrices for the diagrams in each dimension
  # dims are the homological dimensions of diagrams to consider
  # p is a number >=1
  # q is a finite number >= 1
  # distance is the distance metric to use, either "wasserstein" or "fisher"
  # sigma is the positive bandwith for persistence fisher distance

  # create combination of all pairs of diagram group elements and their group indices
  combinations <- do.call(rbind,lapply(X = 1:length(diagram_groups),FUN = function(X){

    distance_pairs <- as.data.frame(t(as.data.frame(utils::combn(x = length(diagram_groups[[X]]),m = 2,simplify = F))))
    distance_pairs$group <- X
    rownames(distance_pairs) <- NULL
    return(distance_pairs[,c(3,1,2)])

  }))

  # determine maximum available number of cores
  num_workers <- parallelly::availableCores(omit = 1)

  # initialize a cluster cl for computing distances between diagrams in parallel
  cl <- parallel::makeCluster(num_workers)
  doParallel::registerDoParallel(cl)

  # export necessary libraries and variables to cl
  parallel::clusterEvalQ(cl,c(library(clue),library(rdist)))
  parallel::clusterExport(cl,c("check_diagram","diagram_distance","diagram_groups","dist_mats","dims","combinations","p"),envir = environment())

  # initialize return vector of statistics, one for each dimension
  statistics <- c()

  # compute loss function and update distance matrices in each dimension dim
  for(dim in dims)
  {

    parallel::clusterExport(cl,"dim")

    d_tots <- foreach::`%dopar%`(obj = foreach::foreach(comb = 1:nrow(combinations),.combine = c),ex = {

      # get group and diagram indices from combinations
      g <- as.numeric(combinations[comb,1])
      d1 <- as.numeric(combinations[comb,2])
      d2 <- as.numeric(combinations[comb,3])

      # get index of dim in dims
      dim_ind <- min(which(dims == dim))

      # if the distance between these two diagrams has not already been computed, compute their distance
      if(dist_mats[dim_ind][[1]][diagram_groups[[g]][[d1]]$ind,diagram_groups[[g]][[d2]]$ind] == -1)
      {
        return(diagram_distance(D1 = diagram_groups[[g]][[d1]]$diag,D2 = diagram_groups[[g]][[d2]]$diag,p = p,dim = dim,distance = distance,sigma = sigma)^q)
      }

      # else return the already stored distance value
      return(dist_mats[dim_ind][[1]][diagram_groups[[g]][[d1]]$ind,diagram_groups[[g]][[d2]]$ind])

    })

    # update the upper triangle of dist_mat to account for new distances just calculated
    for(comb in 1:nrow(combinations))
    {
      # get group and diagram indices
      g <- as.numeric(combinations[comb,1])
      d1 <- as.numeric(combinations[comb,2])
      d2 <- as.numeric(combinations[comb,3])

      # get index of dim in dims
      dim_ind <- min(which(dims == dim))

      # update the correct entry of the current dimension distance matrix with the newly calculated distances
      dist_mats[dim_ind][[1]][diagram_groups[[g]][[d1]]$ind,diagram_groups[[g]][[d2]]$ind] = d_tots[[comb]]

    }

    # append calculated loss statistic to statistics vector
    statistics <- c(statistics,sum(unlist(lapply(X = 1:length(diagram_groups),FUN = function(X){

      # loss statistic is the sum of all within group distances of unique diagram pairs, normalized
      # by the number of those pairs in the group
      return(sum(d_tots[which(combinations$group == X)])/(length(diagram_groups[[X]])*(length(diagram_groups[[X]]) - 1)))

    }))))

  }

  # cleanup
  parallel::stopCluster(cl)

  # return the test statistics and distance matrices in all dimensions
  return(list(statistics = statistics,dist_mats = dist_mats))

}
