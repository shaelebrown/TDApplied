#### DIAGRAM DISTANCE METRICS ####
#' Calculate distance between a pair of persistence diagrams.
#'
#' Calculates the distance between a pair of persistence diagrams, either the output from a \code{\link{diagram_to_df}} function call
#' or from a persistent homology calculation like \code{\link[TDA]{ripsDiag}}/\code{\link[TDAstats]{calculate_homology}}/\code{\link{PyH}},
#' in a particular homological dimension.
#'
#' The most common distance calculations between persistence diagrams
#' are the wasserstein and bottleneck distances, both of which "match" points between
#' their two input diagrams and compute the "loss" of the optimal matching 
#' (see \url{http://www.geometrie.tugraz.at/kerber/kerber_papers/kmn-ghtcpd_journal.pdf} for details). Another 
#' method for computing distances, the Fisher information metric, 
#' converts the two diagrams into distributions
#' defined on the plane, and calculates a distance between the resulting two distributions
#' (\url{https://proceedings.neurips.cc/paper/2018/file/959ab9a0695c467e7caf75431a872e5c-Paper.pdf}).
#' If the `distance` parameter is "fisher" then `sigma` must not be NULL. As noted in the Persistence Fisher paper,
#' there is a fast speed-up approximation which has been implemented from \url{https://github.com/vmorariu/figtree} 
#' and can be accessed by setting the `rho` parameter. Smaller
#' values of `rho` will result in tighter approximations at the expense of longer runtime, and vice versa.
#'
#' @param D1 the first persistence diagram.
#' @param D2 the second persistence diagram.
#' @param dim the non-negative integer homological dimension in which the distance is to be computed, default 0.
#' @param p  a number representing the wasserstein power parameter, at least 1 and default 2.
#' @param distance a string which determines which type of distance calculation to carry out, either "wasserstein" (default) or "fisher".
#' @param sigma either NULL (default) or a positive number representing the bandwidth for the Fisher information metric.
#' @param rho either NULL (default) or a positive number. If NULL then the exact calculation of the Fisher information metric is returned and otherwise a fast approximation, see details.
#'
#' @return the numeric value of the distance calculation.
#' @importFrom rdist cdist
#' @importFrom clue solve_LSAP
#' @import Rcpp
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @references
#' Kerber M, Morozov D and Nigmetov A (2017). "Geometry Helps to Compare Persistence Diagrams." \url{http://www.geometrie.tugraz.at/kerber/kerber_papers/kmn-ghtcpd_journal.pdf}.
#' 
#' Le T, Yamada M (2018). "Persistence fisher kernel: a riemannian manifold kernel for persistence diagrams." \url{https://proceedings.neurips.cc/paper/2018/file/959ab9a0695c467e7caf75431a872e5c-Paper.pdf}.
#' 
#' Vlad I. Morariu, Balaji Vasan Srinivasan, Vikas C. Raykar, Ramani Duraiswami, and Larry S. Davis. Automatic online tuning for fast Gaussian summation. Advances in Neural Information Processing Systems (NIPS), 2008.
#' 
#' @examples
#'
#' if(require("TDA"))
#' {
#'   # create two diagrams
#'   D1 <- TDA::ripsDiag(TDA::circleUnif(n = 20,r = 1),
#'                       maxdimension = 1,maxscale = 2)
#'   D2 <- TDA::ripsDiag(TDA::sphereUnif(n = 20,d = 2,r = 1),
#'                       maxdimension = 1,maxscale = 2)
#'
#'   # calculate 2-wasserstein distance between D1 and D2 in dimension 1
#'   diagram_distance(D1,D2,dim = 1,p = 2,distance = "wasserstein")
#' 
#'   # calculate bottleneck distance between D1 and D2 in dimension 0
#'   diagram_distance(D1,D2,dim = 0,p = Inf,distance = "wasserstein")
#' 
#'   # Fisher information metric calculation between D1 and D2 for sigma = 1 in dimension 1
#'   diagram_distance(D1,D2,dim = 1,distance = "fisher",sigma = 1)
#'   
#'   # repeat but with fast approximation
#'   diagram_distance(D1,D2,dim = 1,distance = "fisher",sigma = 1,rho = 0.001)
#' }

diagram_distance <- function(D1,D2,dim = 0,p = 2,distance = "wasserstein",sigma = NULL,rho = NULL){

  # function to compute the wasserstein/bottleneck/Fisher information metric between two diagrams

  # for standalone usage force D1 and D2 to be data frames if they are the output of a homology calculation
  D1 <- check_diagram(D1,ret = T)
  D2 <- check_diagram(D2,ret = T)
  
  # error check other parameters
  check_param("dim",dim,whole_numbers = T,positive = F,numeric = T,multiple = F,finite = T,non_negative = T)
  check_param("distance",distance)
  check_param("p",p,at_least_one = T,finite = F,non_negative = F,numeric = T,multiple = F)
  if(distance == "fisher")
  {
    check_param("sigma",sigma,positive = T,non_negative = F,numeric = T,finite = T,multiple = F)
    if(!is.null(rho))
    {
      check_param("rho",rho,positive = T,non_negative = T,numeric = T,finite = T,multiple = F)
    }
  }
  
  # subset both diagrams by dimension dim and for birth and death columns
  D1_subset <- D1[which(D1$dimension == dim),1:3]
  D2_subset <- D2[which(D2$dimension == dim),1:3]
  D1_subset <- D1_subset[,2:3]
  D2_subset <- D2_subset[,2:3]
  
  # check if there are any Inf values
  if(length(which(is.infinite(D1_subset[,2]))) > 0 | length(which(is.infinite(D2_subset[,2]))))
  {
    stop("Diagrams must not contain Inf values. This can occur in some persistent homology calculations, like \'alphaComplexDiag\'. Try converting such diagrams to dataframes using the \'diagram_to_df\' function and replacing Inf values with a suitable number (at least the largest finite death value).")
  }

  # create empty diagonals for the persistence diagrams
  diag1 <- D1_subset[0,]
  diag2 <- D2_subset[0,]
  
  # remove diagonal entries from D1_subset and D2_subset
  if(nrow(D1_subset) > 0)
  {
    D1_subset <- D1_subset[which(D1_subset[,1] != D1_subset[,2]),]
  }
  if(nrow(D2_subset) > 0)
  {
    D2_subset <- D2_subset[which(D2_subset[,1] != D2_subset[,2]),]
  }

  # if both subsets are empty then set their distance to 0
  if(nrow(D1_subset) == 0 & nrow(D2_subset) == 0)
  {
    return(0)
  }
  
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
  
  # check if either subset is empty, if so return the norm of the other diagram
  if(distance == "wasserstein")
  {
    if(nrow(diag1) == 0)
    {
      if(is.infinite(p))
      {
        # bottleneck norm
        return(max(unlist(lapply(X = 1:nrow(D2_subset),FUN = function(X){
          
          coord <- (D2_subset[X,1L] + D2_subset[X,2L])/2
          return(max(c(D2_subset[X,2L] - coord,coord - D2_subset[X,1L])))
          
        }))))
      }else
      {
        # wasserstein norm
        return((sum(unlist(lapply(X = 1:nrow(D2_subset),FUN = function(X){
          
          coord <- (D2_subset[X,1L] + D2_subset[X,2L])/2
          return(max(c(D2_subset[X,2L] - coord,coord - D2_subset[X,1L]))^p)
          
        }))))^(1/p))
      }
    }
    if(nrow(diag2) == 0)
    {
      if(is.infinite(p))
      {
        # bottleneck norm
        return(max(unlist(lapply(X = 1:nrow(D1_subset),FUN = function(X){
          
          coord <- (D1_subset[X,1L] + D1_subset[X,2L])/2
          return(max(c(D1_subset[X,2L] - coord,coord - D1_subset[X,1L])))
          
        }))))
      }else
      {
        # wasserstein norm
        return((sum(unlist(lapply(X = 1:nrow(D1_subset),FUN = function(X){
          
          coord <- (D1_subset[X,1L] + D1_subset[X,2L])/2
          return(max(c(D1_subset[X,2L] - coord,coord - D1_subset[X,1L]))^p)
          
        }))))^(1/p))
      }
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
    dist_mat[(nrow(D1_subset) - nrow(diag2) + 1):nrow(D1_subset),(nrow(D2_subset) - nrow(diag1) + 1):nrow(D2_subset)] <- 0
    if(is.finite(p))
    {
      dist_mat <- dist_mat^p
    }
    
    # use the Hungarian algorithm from the clue package to find the minimal weight matching
    best_match <- as.numeric(clue::solve_LSAP(x = dist_mat,maximum = F))
    seq_match <- 1:length(best_match)
    
    # subset best match by removing all pairs between diagonal points
    indices <- cbind(seq_match,best_match)
    indices <- indices[which(indices[,1] <= (nrow(D1_subset) - nrow(diag2)) | indices[,2] <= (nrow(D2_subset) - nrow(diag1))),]
    
    # if p is finite, exponentiate each matched distance and take the p-th root of their sum
    if(is.finite(p))
    {
      return(sum(dist_mat[indices])^(1/p))
    }
    
    # otherwise, return the regular bottleneck distance
    return(max(dist_mat[indices]))
    
  }else
  {
    # compute the persistence Fisher distance
    
    # get all unique points in both diagrams and their diagonal projections
    theta <- rbind(D1_subset,D2_subset)
    
    # exact calculation, quadratic runtime
    if(is.null(rho))
    {
      rho_1 <- unlist(lapply(X = 1:nrow(theta),FUN = function(X){
        
        x <- as.numeric(theta[X,])
        sum <- 0
        for(i in 1:nrow(D1_subset))
        {
          # compute mvn normal pdf values and sum up
          u <- as.numeric(D1_subset[i,])
          sum <- sum + exp(-1*((x[[1]]-u[[1]])^2 + (x[[2]]-u[[2]])^2)/(2*sigma^2))/(sqrt(2*pi)*sigma)
        }
        return(sum)
        
      }))
      rho_2 <- unlist(lapply(X = 1:nrow(theta),FUN = function(X){
        
        x <- as.numeric(theta[X,])
        sum <- 0
        for(i in 1:nrow(D2_subset))
        {
          # compute mvn normal pdf values and sum up
          u <- as.numeric(D2_subset[i,])
          sum <- sum + exp(-1*((x[[1]]-u[[1]])^2 + (x[[2]]-u[[2]])^2)/(2*sigma^2))/(sqrt(2*pi)*sigma)
        }
        return(sum)
        
      }))
    }else
    {
      # approximation, linear runtime
      # must convert tables into vectors of length 2*nrow
      # of the form c(first_row_birth,first_row_death,second_row_birth,second_row_death,...)
      D1_subset <- as.matrix(D1_subset)
      D2_subset <- as.matrix(D2_subset)
      theta <- as.matrix(theta)
      
      D1_subset <- as.vector(t(D1_subset))
      D2_subset <- as.vector(t(D2_subset))
      theta <- as.vector(t(theta))
      
      # this is the results vector which will be
      # modified and returned
      G <- rep(0,length(theta)/2)
      
      # use figtree c++ code to speed up calculations
      rho_1 <- as.numeric(figtree(X = D1_subset,h = sqrt(2)*sigma,Q = rep(1,length(D1_subset)/2)/(length(D1_subset)/2),Y = theta,epsilon = rho,G = G))
      rho_2 <- as.numeric(figtree(X = D2_subset,h = sqrt(2)*sigma,Q = rep(1,length(D2_subset)/2)/(length(D2_subset)/2),Y = theta,epsilon = rho,G = G))
    }
    
    if(length(which(rho_1 == rho_2)) == length(rho_1)) # same diagrams
    {
      return(0)
    }
    
    rho_1 <- rho_1/sum(rho_1)
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
#' Compute a distance matrix from a list of persistence diagrams.
#'
#' Calculate the distance matrix \eqn{d} for either a single list of persistence diagrams \eqn{(D_1,D_2,\dots,D_n)}, i.e. \eqn{d[i,j] = d(D_i,D_j)}, 
#' or between two lists, \eqn{(D_1,D_2,\dots,D_n)} and \eqn{(D'_1,D'_2,\dots,D'_n)}, \eqn{d[i,j] = d(D_i,D'_j)}, in parallel.
#'
#' Distance matrices of persistence diagrams are used in downstream analyses, like in the 
#' \code{\link{diagram_mds}}, \code{\link{permutation_test}} and \code{\link{diagram_ksvm}} functions. 
#' If `distance` is "fisher" then `sigma` must not be NULL. Since the matrix is computed sequentially when
#' approximating the Fisher information metric this is only recommended when the persistence diagrams
#' contain many points and when the number of available cores is small.
#'
#' @param diagrams a list of persistence diagrams, either the output of persistent homology calculations like \code{\link[TDA]{ripsDiag}}/\code{\link[TDAstats]{calculate_homology}}/\code{\link{PyH}}, or \code{\link{diagram_to_df}}.
#' @param other_diagrams either NULL (default) or another list of persistence diagrams to compute a cross-distance matrix.
#' @param dim the non-negative integer homological dimension in which the distance is to be computed, default 0.
#' @param distance a character determining which metric to use, either "wasserstein" (default) or "fisher".
#' @param p a number representing the wasserstein power parameter, at least 1 and default 2.
#' @param sigma a positive number representing the bandwidth of the Fisher information metric, default NULL.
#' @param rho an optional positive number representing the heuristic for Fisher information metric approximation, see \code{\link{diagram_distance}}. Default NULL. If not NULL then matrix is calculated sequentially, but functions in the "exec" directory
#'            of the package can be loaded to calculate distance matrices in parallel with approximation.
#' @param num_workers the number of cores used for parallel computation, default is one less than the number of cores on the machine.
#'
#' @return the numeric distance matrix.
#' @export
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom parallel makeCluster stopCluster clusterExport
#' @importFrom parallelly availableCores
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom iterators iter
#' @examples
#'
#' if(require("TDA") & require("TDAstats"))
#' {
#'   # create two diagrams
#'   D1 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                      dim = 0,threshold = 2)
#'   D2 <- TDAstats::calculate_homology(TDA::circleUnif(n = 10,r = 1),
#'                                      dim = 0,threshold = 2)
#'   g <- list(D1,D2)
#'
#'   # calculate their distance matrix in dimension 0 with the persistence Fisher metric
#'   # using 2 cores
#'   D <- distance_matrix(diagrams = g,dim = 0,distance = "fisher",sigma = 1,num_workers = 2)
#'   
#'   # calculate their distance matrix in dimension 0 with the approximate persistence Fisher metric
#'   # using 2 cores
#'   D <- distance_matrix(diagrams = g,dim = 0,distance = "fisher",sigma = 1,rho = 0.001,
#'                        num_workers = 2)
#'
#'   # calculate their distance matrix in dimension 0 with the 2-wasserstein metric 
#'   # using 2 cores
#'   D <- distance_matrix(diagrams = g,dim = 0,distance = "wasserstein",p = 2,num_workers = 2)
#' 
#'   # now do the cross distance matrix, which is the same as the previous
#'   D_cross <- distance_matrix(diagrams = g,other_diagrams = g,
#'                              dim = 0,distance = "wasserstein",
#'                              p = 2,num_workers = 2)
#' }

distance_matrix <- function(diagrams,other_diagrams = NULL,dim = 0,distance = "wasserstein",p = 2,sigma = NULL,rho = NULL,num_workers = parallelly::availableCores(omit = 1)){
  
  # calculate a (cross) distance matrix of persistence diagrams in parallel (if desired)
  
  # set internal variables to NULL to avoid build issues
  r <- NULL
  X <- NULL
  
  # error check diagrams argument
  check_param("diagrams",diagrams,min_length = 1)
  diagrams <- all_diagrams(diagram_groups = list(diagrams),inference = "independence")[[1]]
  if(!is.null(other_diagrams))
  {
    check_param("other_diagrams",other_diagrams,min_length = 1)
    other_diagrams <- all_diagrams(diagram_groups = list(other_diagrams),inference = "independence")[[1]]
  }
  
  # check other parameters
  check_param("p",p,finite = F,at_least_one = T,numeric = T,multiple = F)
  check_param("dim",dim,whole_numbers = T,non_negative = T,positive = F,numeric = T,multiple = F)
  check_param("distance",distance)
  if(distance == "fisher")
  {
    check_param("sigma",sigma,positive = T,numeric = T,finite = T,multiple = F)
    if(!is.null(rho))
    {
      check_param("rho",rho,positive = T,non_negative = T,numeric = T,finite = T,multiple = F)
    }
  }
  
  # error check num_workers argument
  check_param("num_workers",num_workers,whole_numbers = T,at_least_one = T,finite = T,numeric = T,multiple = F)
  if(num_workers > parallelly::availableCores())
  {
    warning("num_workers is greater than the number of available cores - setting to maximum value less one.")
    num_workers <- parallelly::availableCores(omit = 1)
  }

  # set up cluster
  m = length(diagrams)
  cl <- parallel::makeCluster(num_workers)
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl,varlist = c("diagram_distance","check_diagram","check_param","figtree"),envir = environment())
  
  # if approximation to Fisher information metric is used, run sequentially
  # otherwise in parallel
  if(distance == "fisher" & !is.null(rho))
  {
    foreach_func <- foreach::`%do%`
  }else
  {
    foreach_func <- foreach::`%dopar%`
  }
  
  # catch errors during parallel calculations to make sure
  # clusters are closed
  tryCatch(expr = {
    
    if(is.null(other_diagrams))
    {
      
      # not cross distance matrix, only need to compute the upper diagonal
      # since the matrix is symmetric
      d <- matrix(data = 0,nrow = length(diagrams),ncol = length(diagrams))
      u <- which(upper.tri(d),arr.ind = T)
      R <- lapply(X = 1:nrow(u),FUN = function(X){
        
        return(list(diagrams[[u[[X,1]]]],diagrams[[u[[X,2]]]]))
        
      })
      
      # remove diagrams to preserve memory
      rm(diagrams)
      
      # calculate distances in parallel, export TDApplied to nodes
      d_off_diag <- foreach_func(obj = foreach::foreach(r = R,.combine = c,.packages = c("clue","rdist")),ex = {diagram_distance(D1 = r[[1]],D2 = r[[2]],dim = dim,distance = distance,p = p,sigma = sigma,rho = rho)})
      
      # store results in matrix
      d[upper.tri(d)] <- d_off_diag
      d[which(upper.tri(d),arr.ind = T)[,c("col","row")]] <- d_off_diag
      diag(d) <- rep(0,nrow(d))
      
    }else
    {

      # cross distance matrix, need to compute all entries
      d <- matrix(data = 0,nrow = length(diagrams),ncol = length(other_diagrams))
      u <- expand.grid(1:length(diagrams),1:length(other_diagrams))
      R <- lapply(X = 1:nrow(u),FUN = function(X){
        
        return(list(diagrams[[u[X,1]]],other_diagrams[[u[X,2]]]))
        
      })
      
      # remove diagrams and other_diagrams to preserve memory
      rm(list = c("diagrams","other_diagrams"))
      
      # store distance calculations in matrix
      d[as.matrix(u)] <- foreach_func(foreach::foreach(r = R,.combine = cbind,.export = c("diagram_distance","check_diagram","check_param","figtree"),.packages = c("clue","rdist")),ex = {diagram_distance(D1 = r[[1]],D2 = r[[2]],dim = dim,distance = distance,p = p,sigma = sigma,rho = rho)})
      
    }
    
  }, warning = function(w){warning(w)},
  error = function(e){stop(e)},
  finally = {
    # close cluster
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
    
  })
  
  return(d)
  
}

#### LOSS FUNCTION FOR GROUPS OF PERSISTENCE DIAGRAMS ####
#' Turner loss function for a list of groups (lists) of persistence diagrams.
#'
#' An internal function to calculate the normalized sum of within-group exponentiated distances 
#' between pairs of persistence diagrams (stored as data frames)
#' for an arbitrary number of groups in parallel. Note that this function may run
#' into memory issues for large numbers of diagrams.
#' 
#' The Turner loss function is described in Robinson and Turner 2017
#' (\url{https://link.springer.com/article/10.1007/s41468-017-0008-7}), and is used
#' in the `permutation_test` function to describe how well-separated a particular
#' grouping of persistence diagrams is. When the `distance` parameter is "fisher",
#' `sigma` must not be NULL.
#'
#' @param diagram_groups groups (lists/vectors) of persistence diagrams, stored as lists of a data frame and
#'                          an index of the diagram in all the diagrams across all groups.
#' @param dist_mats distance matrices between all possible pairs of persistence diagrams across and within groups
#'                      storing the current distances which have been pre-computed.
#' @param dims a numeric vector of which homological dimensions in which the loss function is to be computed.
#' @param p a number representing the wasserstein parameter, at least 1, and if Inf then the bottleneck distance is calculated.
#' @param q a finite number at least 1.
#' @param distance a string which determines which type of distance calculation to carry out, either "wasserstein" (default) or "fisher".
#' @param sigma the positive bandwidth for the persistence Fisher distance.
#' @param rho the approximation heuristic for Fisher information metric, results in sequential computation.
#' @param num_workers the number of cores used for parallel computation.
#' @param group_sizes for when using precomputed distance matrices.
#'
#' @importFrom parallel makeCluster clusterExport stopCluster
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom utils combn
#' @author Shael Brown - \email{shaelebrown@@gmail.com}
#' @keywords internal
#' @references
#' Robinson T, Turner K (2017). "Hypothesis testing for topological data analysis." \url{https://link.springer.com/article/10.1007/s41468-017-0008-7}.
#' @return the numeric value of the Turner loss function.

loss <- function(diagram_groups,dist_mats,dims,p,q,distance,sigma,rho,num_workers,group_sizes){

  # function to compute the F_{p,q} loss between groups of diagrams
  # diagram_groups are the (possibly permuted) groups of diagrams
  
  # if approximation to Fisher information metric is used, run sequentially
  # otherwise in parallel
  if(distance == "fisher" & !is.null(rho))
  {
    foreach_func <- foreach::`%do%`
  }else
  {
    foreach_func <- foreach::`%dopar%`
  }

  if(is.numeric(diagram_groups[[1]][[1]]) == F)
  {
    # distance matrices were not precomputed
    
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
    
    # export necessary functions to cl
    parallel::clusterExport(cl,c("check_diagram","diagram_distance","check_param","figtree"),envir = environment())
    
    tryCatch(expr = {
      
      # initialize return vector of statistics, one for each dimension
      statistics <- c()
      
      # compute loss function and update distance matrices in each dimension dim
      for(dim in dims)
      {
        
        d_tots <- foreach_func(obj = foreach::foreach(comb = 1:nrow(combinations),.combine = c,.packages = c("clue","rdist")),ex = {
          
          # get group and diagram indices from combinations
          g <- as.numeric(combinations[comb,1])
          d1 <- as.numeric(combinations[comb,2])
          d2 <- as.numeric(combinations[comb,3])
          
          # get index of dim in dims
          dim_ind <- min(which(dims == dim))
          
          # if the distance between these two diagrams has not already been computed, compute their distance
          if(dist_mats[dim_ind][[1]][diagram_groups[[g]][[d1]]$ind,diagram_groups[[g]][[d2]]$ind] == -1)
          {
            res <- diagram_distance(D1 = diagram_groups[[g]][[d1]]$diag,D2 = diagram_groups[[g]][[d2]]$diag,p = p,dim = dim,distance = distance,rho = rho,sigma = sigma)^q
          }else
          {
            # else return the already stored distance value
            res <- dist_mats[dim_ind][[1]][diagram_groups[[g]][[d1]]$ind,diagram_groups[[g]][[d2]]$ind]
          }
          res
          
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
      
    },
    warning = function(w){warning(w)},
    error = function(e){stop(e)},
    finally = {
      
      # cleanup
      doParallel::stopImplicitCluster()
      parallel::stopCluster(cl)
      
    }) 
  }else
  {
    # dist_mats are supplied in full
    statistics <- unlist(lapply(X = dims,FUN = function(X){
      
      # sum starting at 0
      s <- 0
      
      # add to sum for each group
      for(g in 1:length(group_sizes))
      {
        # get all pairs of indices in the group g
        inds <- diagram_groups[[g]]
        inds <- expand.grid(inds,inds)
        inds <- inds[which(as.numeric(inds[,1]) < as.numeric(inds[,2])),]
        colnames(inds) <- NULL
        rownames(inds) <- NULL
        inds <- as.matrix(inds)
        
        # add normalized sum of within group distances
        s <- s + (sum((dist_mats[[which(dims == X)[[1]]]][unlist(inds)])^q))/(group_sizes[[g]]*(group_sizes[[g]] - 1))
      }
      return(s)
      
    }))
  }

  # return the test statistics and distance matrices in all dimensions
  return(list(statistics = statistics,dist_mats = dist_mats))

}
