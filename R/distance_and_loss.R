#### DIAGRAM DISTANCE METRICS ####
#' Calculate distances between pairs of persistence diagrams
#'
#' Calculates the distance between a pair of persistce diagrams, stored as
#' data frames (as the output from diagram_to_df)
#' in a particular homological dimension. Different TDA sources define distances
#' differently, and this function has functionality to compute distances like
#' in the R package TDA (based on the C++ library Dionysus, see
#' <https://mrzv.org/software/dionysus2/>) or like in the
#' original paper for inference of persistence diagrams by Robinson and Turner in 2017
#' <https://link.springer.com/article/10.1007/s41468-017-0008-7>.
#'
#' The `D1` and `D2` parameters should be persistence diagrams, outputted
#' from a homology calculation in the package TDA, or such a
#' persistence diagram converted to a data frame via the function diagram_to_df.
#' The `dim` parameter should be a positive finite integer.
#' The `p` parameter should be a positive integer or Inf. The `distance` parameter
#' should be a string, either "wasserstein" or "Turner".
#'
#' @param D1 first persistence diagram, either computed from TDA or converted to a data frame
#' @param D2 second persistence diagram, either computed from TDA or converted to a data frame
#' @param dim homological dimension in which the distance is to be computed
#' @param p  matching distance parameter
#' @param distance string which determines which type of distance calculation to carry out
#'
#' @return numeric value of the distance calculation
#' @importFrom rdist cdist
#' @importFrom clue solve_LSAP
#' @export
#' @examples
#'
#' # create two diagrams with package TDA based on 2D Gaussians
#' diag1 <- TDA::ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),y = rnorm(100,mean = 0,sd = 1)),maxscale = 1,maxdimension = 1)
#' diag2 <- TDA::ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),y = rnorm(100,mean = 0,sd = 1)),maxscale = 1,maxdimension = 1)
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

diagram_distance <- function(D1,D2,dim,p,distance){

  # function to compute the wasserstein/bottleneck metric between two diagrams
  # D1 and D2 are diagrams stored as data frames
  # dim is the dimension to subset
  # p is the power of the wasserstein distance, p >= 1
  # distance is either "wasserstein" (default) or "Turner"
  # if p == Inf or distance == "wasserstein" then the underlying matching distance is
  # the bottleneck distance

  # for standalone usage force D1 and D2 to be data frames if they are the output of a homology calculation
  if(length(class(D1)) != 1 && class(D1) != "data.frame")
  {
    if(is.list(D1) && class(D1[[1]]) == "diagram")
    {
      # D1 is the output from a TDA calculation
      D1 <- diagram_to_df(D1)
    }else
    {
      stop("D1 must be the output of either a TDA computation.")
    }

    # error check D1
    check_diagram(D1)

  }

  if(length(class(D2)) != 1 && class(D2) != "data.frame")
  {
    if(is.list(D2) && class(D2[[1]]) == "diagram")
    {
      # D2 is the output from a TDA calculation
      D2 <- TDA_diagram_to_df(D2)
    }else
    {
      stop("D2 must be the output of either a TDA computation.")
    }

    # error check for D2
    check_diagram(D2)

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

  # compute the distance matrix between rows of D1_subset and D2_subset
  if(is.finite(p) & distance == "Turner")
  {
    # compute the p-wasserstein distance for matching pairs
    dist_mat <- as.matrix(rdist::cdist(D1_subset,D2_subset,metric = "minkowski",p = p))
  }else
  {
    # compute the bottleneck distance for matching pairs
    dist_mat <- as.matrix(rdist::cdist(D1_subset,D2_subset,metric = "maximum"))
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
    return((sum(dist_mat[indices]^(p)))^(1/p))
  }

  # otherwise, return the regular bottleneck distance
  return(max(dist_mat[indices]))

}

#### LOSS FUNCTION FOR GROUPS OF PERSISTENCE DIAGRAMS ####
#' Calculate Turner loss function for groups of persistence diagrams
#'
#' Calculates the pairs of distances between all persistence diagrams (stored as data frames)
#' with an arbitrary number of groups. The loss function is described in Robinson and Turner in 2017
#' <https://link.springer.com/article/10.1007/s41468-017-0008-7>, but mathematically we
#' compute the distances of each within-group pair of persistence diagrams and exponentiate
#' each distance by q and take the q-th root of the sum.
#'
#' The `diagram_groups` parameter should be a list or vector of persistence diagrams stored as data frames.
#' The `dist_mats` parameter should be a list, with one element for each element in the parameter `dims`,
#' which stores a matrix of distance calculations (and -1 for distance calculations yet to be completed).
#' The `dims` parameter is a vector of non-negative whole numbers, representing the homological dimensions
#' to calculate the loss function in. The `p` parameter should be a number at least 1 and possibly Inf.
#' The `q` parameter should be a finite number at least 1. The `distance` parameter should be a string
#' either "wasserstein" or "Turner".
#'
#' @param diagram_groups groups (lists/vectors) of persistence diagrams, stored as lists of a data frame and
#'                          an index of the diagram in all the diagrams across all groups
#' @param dist_mats distance matrices between all possible pairs of persistence diagrams across and within groups
#'                      storing the current distances which have been precomputed.
#' @param dims which homological dimensions in which the loss function is to be computed
#' @param p positive wasserstein parameter, if Inf then the bottleneck distance
#' @param q finite exponent at least 1
#' @param distance string which determines which type of distance calculation to carry out
#'
#' @importFrom parallel makeCluster detectCores clusterEvalQ clusterExport stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom utils combn
#' @return numeric value of the loss function
#' @examples
#'
#' # create two groups of diagrams, based on 2D Gaussians, with package TDA
#' g1 <- lapply(X = 1:3,FUN = function(X){
#'
#' diag <- TDA::ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),y = rnorm(100,mean = 0,sd = 1)),maxscale = 1,maxdimension = 1)
#' df <- diagram_to_df(d = diag)
#' return(list(diag = df,ind = X))
#' })
#'
#' g2 <- lapply(X = 1:3,FUN = function(X){
#'
#' diag <- TDA::ripsDiag(data.frame(x = rnorm(100,mean = 0,sd = 1),y = rnorm(100,mean = 0,sd = 1)),maxscale = 1,maxdimension = 1)
#' df <- diagram_to_df(d = diag)
#' return(list(diag = df,ind = X + 3))
#' })
#'
#'
#' # compute Turner loss function with p,q = 2 in dimension 1 with other parameters set as defaults
#' example_loss <- loss(diagram_groups = list(g1,g2),dist_mats = list(matrix(data = -1,nrow = 6,ncol = 6)),p = 2,q = 2,distance = "Turner")

loss <- function(diagram_groups,dist_mats,dims,p,q,distance){

  # function to compute the F_{p,q} loss between groups of diagrams
  # diagram_groups are the (possibly permuted) groups of diagrams
  # dist mats stores the current distance matrices for the diagrams in each dimension
  # dims are the homological dimensions of diagrams to consider
  # p is a number >=1
  # q is a finite number >= 1
  # distance is the distance metric to use, either "wasserstein" (default) or "Turner"

  # create combination of all pairs of diagram group elements and their group indices
  combinations <- do.call(rbind,lapply(X = 1:length(diagram_groups),FUN = function(X){

    distance_pairs <- as.data.frame(t(as.data.frame(utils::combn(x = length(diagram_groups[[X]]),m = 2,simplify = F))))
    distance_pairs$group <- X
    rownames(distance_pairs) <- NULL
    return(distance_pairs[,c(3,1,2)])

  }))

  # initialize a cluster cl for computing distances between diagrams in parallel
  cl <- parallel::makeCluster(parallel::detectCores() - 1)
  doParallel::registerDoParallel(cl)

  # export necessary libraries and variables to cl
  parallel::clusterEvalQ(cl,c(library(clue),library(rdist)))
  parallel::clusterExport(cl,c("diagram_distance","diagram_groups","dist_mats","dims","combinations","p"),envir = environment())

  # initialize return vector of statistics, one for each dimension
  statistics <- c()

  # compute loss function and update distance matrices in each dimension dim
  for(dim in dims)
  {

    parallel::clusterExport(cl,"dim",envir = environment())

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
        return(diagram_distance(D1 = diagram_groups[[g]][[d1]]$diag,D2 = diagram_groups[[g]][[d2]]$diag,p = p,dim = dim,distance = distance)^q)
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
