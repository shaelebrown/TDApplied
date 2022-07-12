#### DIAGRAM DISTANCE METRICS ####
#' Calculate distances between pairs of persistence diagrams
#'
#' Calculates the distance between a pair of persistence diagrams in a particular homological dimension,
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
#' # load three diagrams
#' D1 <- generate_TDAML_test_data(1,0,0)
#' D2 <- generate_TDAML_test_data(0,1,0)
#' D3 <- generate_TDAML_test_data(0,0,1)
#'
#' # calculate 2-wasserstein distance between D1 and D2
#' diagram_distance(D1,D2,dim = 0,p = 2,distance = "wasserstein")
#' 
#' # calculate bottleneck distance between D1 and D3
#' diagram_distance(D1,D3,dim = 0,p = Inf,distance = "wasserstein")
#' 
#' # Fisher information metric calculation between D1 and D2 for sigma = 1
#' diagram_distance(D1,D2,dim = 0,distance = "fisher",sigma = 1)
#'
#' # Fisher information metric calculation between D1 and D3 for sigma = 2
#' diagram_distance(D1,D3,dim = 0,distance = "fisher",sigma = 2)

diagram_distance <- function(D1,D2,dim,p = 2,distance = "wasserstein",sigma = NULL){

  # function to compute the wasserstein/bottleneck/Fisher information metric between two diagrams
  # D1 and D2 are diagrams, possibly stored as data frames
  # dim is the dimension to subset
  # p is the power of the wasserstein distance, p >= 1, default 2
  # distance is either "wasserstein" (default) or "fisher"
  # sigma is the positive bandwidth for the Fisher information metric, default NULL

  # for standalone usage force D1 and D2 to be data frames if they are the output of a homology calculation
  D1 <- check_diagram(D1,ret = T)
  D2 <- check_diagram(D2,ret = T)
  
  # error check other parameters
  check_param("dim",dim,whole_numbers = T,positive = F)
  check_param("distance",distance)
  check_param("p",p,at_least_one = T,finite = F,non_negative = F)
  if(distance == "fisher")
  {
    check_param("sigma",sigma,positive = T,non_negative = F)
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
    
    # return arc cos of dot product of elementwise square root
    norm <- as.numeric(sqrt(rho_1) %*% sqrt(rho_2))
    if(norm > 1)
    {
      norm <- 1
    }
    if(norm < -1)
    {
      norm <- -1
    }
    return(acos(norm))
    
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
#' `t` is the positive scale parameter for the persistence Fisher kernel. `num_workers` is the
#' number of cores used for parallel computation.
#'
#' @param diagrams the list of persistence diagrams, either the output from TDA calculations or the \code{\link{diagram_to_df}} function.
#' @param other_diagrams either NULL (default) or another list of persistence diagrams to compute a cross-distance matrix.
#' @param dim the homological dimension in which the distance is to be computed.
#' @param distance a character determining which metric to use, either "wasserstein" (default) or "fisher".
#' @param p the positive wasserstein power, default 2.
#' @param sigma a positive number representing the bandwith of the Fisher information metric, default NULL.
#' @param num_workers the number of cores used for parallel computation, default is one less the number of cores on the machine.
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
#' # load three diagrams
#' D1 <- generate_TDAML_test_data(1,0,0)
#' D2 <- generate_TDAML_test_data(0,1,0)
#' D3 <- generate_TDAML_test_data(0,0,1)
#' g <- list(D1,D2,D3)
#'
#' # calculate their distance matrix in dimension 1 with the 2-wasserstein metric using 2 cores
#' D <- distance_matrix(diagrams = g,dim = 1,distance = "wasserstein",p = 2,num_workers = 2)
#' 
#' # now do the cross distance matrix, should be the same as the original
#' D_cross <- distance_matrix(diagrams = g,other_diagrams = g,dim = 1,distance = "wasserstein",p = 2,num_workers = 2)

distance_matrix <- function(diagrams,other_diagrams = NULL,dim = 0,distance = "wasserstein",p = 2,sigma = NULL,num_workers = parallelly::availableCores(omit = 1)){
  
  # set internal variables to NULL to avoid build issues
  r <- NULL
  X <- NULL
  
  # error check diagrams argument
  check_param("diagrams",diagrams,numeric = F,multiple = T)
  diagrams <- all_diagrams(diagram_groups = list(diagrams),inference = "independence")[[1]]
  if(!is.null(other_diagrams))
  {
    check_param("other_diagrams",other_diagrams,numeric = F,multiple = T)
    other_diagrams <- all_diagrams(diagram_groups = list(other_diagrams),inference = "independence")[[1]]
  }
  
  # check other parameters
  check_param("p",p,finite = F,at_least_one = T)
  check_param("dim",dim,whole_numbers = T,non_negative = T,positive = F)
  check_param("distance",distance)
  if(distance == "fisher")
  {
    check_param("sigma",sigma,positive = T)
  }
  
  # error check num_workers argument
  check_param("num_workers",num_workers,whole_numbers = T,at_least_one = T)
  if(num_workers > parallelly::availableCores())
  {
    warning("num_workers is greater than the number of available cores - setting to maximum value.")
    num_workers <- parallelly::availableCores()
  }

  # compute distance matrix in parallel
  m = length(diagrams)
  cl <- parallel::makeCluster(num_workers)
  doParallel::registerDoParallel(cl)
  parallel::clusterEvalQ(cl,c(library(clue),library(rdist)))
  parallel::clusterExport(cl,c("diagram_distance","check_diagram","check_param"),envir = environment())
  force(diagrams) # required for parallel computation in this environment
  
  if(is.null(other_diagrams))
  {
    d <- matrix(data = 0,nrow = m,ncol = m)
    d_off_diag <- foreach::`%dopar%`(obj = foreach::foreach(r = iterators::iter(which(upper.tri(d),arr.ind = T),by = 'row'),.combine = c),ex = {
      
      return(diagram_distance(D1 = diagrams[[r[[1]]]],D2 = diagrams[[r[[2]]]],dim = dim,p = p,distance = distance,sigma = sigma))
      
    })
    d[upper.tri(d)] <- d_off_diag
    d[which(upper.tri(d),arr.ind = T)[,c("col","row")]] <- d_off_diag
    diag(d) <- rep(0,m)
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
#' Fisher distance. `num_workers` is the number of cores used for parallel computation.
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
#' @param num_workers the number of cores used for parallel computation.
#'
#' @importFrom parallel makeCluster clusterEvalQ clusterExport stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom utils combn
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @references
#' Robinson T, Turner K (2017). "Hypothesis testing for topological data analysis." \url{https://link.springer.com/article/10.1007/s41468-017-0008-7}.
#' @return the numeric value of the Turner loss function.

loss <- function(diagram_groups,dist_mats,dims,p,q,distance,sigma,num_workers){

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

  # initialize a cluster cl for computing distances between diagrams in parallel
  cl <- parallel::makeCluster(num_workers)
  doParallel::registerDoParallel(cl)

  # export necessary libraries and variables to cl
  parallel::clusterEvalQ(cl,c(library(clue),library(rdist)))
  parallel::clusterExport(cl,c("check_diagram","diagram_distance","diagram_groups","dist_mats","dims","combinations","p","check_param"),envir = environment())

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
