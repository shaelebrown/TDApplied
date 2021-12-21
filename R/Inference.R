# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


# Script for Statistical Inference using Persistent Hom
# not based on functions in the TDA package, because
# it uses the Dionysus wasserstein distance, which is by default
# an infinity-norm distance

# IMPORT LIBRARIES ----
library(TDA)
library(TDAStats)
library(parallel)
library(doParallel)
library(foreach)
library(clue)
library(rdist)

# FUNCTIONS ----
TDA_diagram_to_df <- function(d){
  # function to convert d to a data frame with standardized colnames names
  # d is a diagram from library TDA

  d = d[[1]]
  class(d) = "matrix"
  d = as.data.frame(d)
  colnames(d) = c("dimension","birth","death")

  return(d)
}

TDAStats_diagram_to_df <- function(d){
  # function to convert d to a dataframe
  # d is a diagram from library TDAStats

  return(as.data.frame(d))
}

check_diagram <- function(d){
  # error checks for a diagram
  # d is a diagram, stored as a dataframe

  if(nrow(d) == 0)
  {
    stop("Every diagram must be non-empty.")
  }

  if(ncol(d) != 3)
  {
    stop("Every diagram must have three columns.")
  }

  if(class(d[,1]) != "numeric" | class(d[,2]) != "numeric" | class(d[,3]) != "numeric")
  {
    stop("Diagrams must have numeric columns.")
  }

  if(!all.equal(d[,1],as.integer(d[,1])))
  {
    stop("Homology dimensions must be whole numbers.")
  }

  if(length(which(d[,1] < 0)) > 0)
  {
    stop("Homology dimensions must be >= 0.")
  }

  if(length(which(d[,2] < 0)) > 0 | length(which(d[,3] < 0)) > 0)
  {
    stop("Birth and death radii must be >= 0.")
  }

  if(length(which(complete.cases(d))) != nrow(d))
  {
    stop("Diagrams can't have missing values.")
  }

}

all_diagrams <- function(diagram_groups,lib){
  # function to make sure all groups are lists or vectors of diagrams, to convert the diagrams to dataframes
  # and to error check each diagram
  # groups is a vector or list of vectors or lists of diagrams
  # lib is either "TDA" or "TDAStats

  # compute cumulative sums of groups lengths in order to correctly compute diagram indices
  csum_group_sizes = cumsum(unlist(lapply(diagram_groups,FUN = length)))
  csum_group_sizes = c(0,csum_group_sizes)

  for(g in 1:length(diagram_groups))
  {
    for(diag in 1:length(diagram_groups[[g]]))
    {
      if(lib == "TDA")
      {
        if(class(diagram_groups[[g]][[diag]][[1]]) != "diagram")
        {
          stop(paste0("Every diagram must be the output from a homology calculation from ",lib,"."))
        }else
        {
          diagram_groups[[g]][[diag]] = list(diag = TDA_diagram_to_df(diagram_groups[[g]][[diag]]),ind = csum_group_sizes[g] + diag)
        }
      }else
      {
        if(class(diagram_groups[[g]][[diag]])[[1]] != "matrix" & class(diagram_groups[[g]][[diag]])[[2]] != "array")
        {
          stop(paste0("Every diagram must be the output from a homology calculation from ",lib,"."))
        }else
        {
          diagram_groups[[g]][[diag]] = list(diag = TDAStats_diagram_to_df(diagram_groups[[g]][[diag]]),ind = csum_group_sizes[g] + diag)
        }
      }

      check_diagram(diagram_groups[[g]][[diag]]$diag)
    }
  }

  return(diagram_groups)

}

check_params <- function(iterations,p,q,dims,paired,lib,distance){
  # error checks on the parameters for the main function

  if(!is.numeric(iterations))
  {
    stop("Iterations must be numeric.")
  }

  if(iterations != as.integer(iterations))
  {
    stop("Iterations must be a whole number.")
  }

  if(iterations <= 1 | is.infinite(iterations))
  {
    stop("Iterations must be a finite number at least 1.")
  }

  if(!as.numeric(p))
  {
    stop("p must be numeric.")
  }

  if(p < 1)
  {
    stop("p must be at least 1.")
  }

  if(!is.numeric(q))
  {
    stop("q must be numeric.")
  }

  if(q < 1 | is.infinite(q))
  {
    stop("q must be a finite number at least 1.")
  }

  if(!is.numeric(dims))
  {
    stop("dims must be a numeric vector or value.")
  }

  if(!all.equal(dims,as.integer(dims)))
  {
    stop("dims must be whole numbers.")
  }

  if(length(which(dims < 0 | is.infinite(dims))) > 0)
  {
    stop("dims must be positive and finite.")
  }

  if(!is.logical(paired))
  {
    stop("paired must be T or F.")
  }

  if(!is.character(lib) | lib %in% c("TDA","TDAStats") == F)
  {
    stop("lib must be a single character, either TDA or TDAStats.")
  }

  if(!is.character(distance) | distance %in% c("wasserstein","Turner") == F)
  {
    stop("distance must be a single character, either wasserstein or Turner.")
  }

}

diagram_distance <- function(D1,D2,dim,p,distance){
  # function to compute the wasserstein/bottleneck metric between two diagrams
  # D1 and D2 are diagrams stored as data frames
  # dim is the dimension to subset
  # p is the power of the wasserstein distance, p >= 1

  # for standalone usage force D1 and D2 to be dataf rames if they are the output of a homology calculation
  if(class(D1) != "data.frame")
  {
    if(is.list(D1) & class(D1[[1]]) == "diagram")
    {
      # D1 is the output from a TDA calculation
      D1 = TDA_diagram_to_df(D1)
    }else
    {
      if(class(D1)[[1]] == "matrix" & class(D1)[[2]] == "array")
      {
        # D1 is the output from a TDAStats calculation
        D1 = TDAStats_diagram_to_df(D1)
      }else
      {
        stop("D1 must be the output of either a TDA or TDAStats computation.")
      }
    }
    check_diagram(D1)
  }

  if(class(D2) != "data.frame")
  {
    if(is.list(D2) & class(D2[[1]]) == "diagram")
    {
      # D2 is the output from a TDA calculation
      D2 = TDA_diagram_to_df(D2)
    }else
    {
      if(class(D2)[[1]] == "matrix" & class(D2)[[2]] == "array")
      {
        # D2 is the output from a TDAStats calculation
        D2 = TDAStats_diagram_to_df(D2)
      }else
      {
        stop("D2 must be the output of either a TDA or TDAStats computation.")
      }
    }
    check_diagram(D2)
  }

  # subset both diagrams by dimension X
  D1_subset = D1[which(D1$dimension == dim),]
  D2_subset = D2[which(D2$dimension == dim),]
  D1_subset = D1_subset[,2:3]
  D2_subset = D2_subset[,2:3]

  # create empty diagonals for the persistence landscapes
  diag1 = D1_subset[0,]
  diag2 = D2_subset[0,]

  # if both subsets are empty then set their distance to 0
  if(nrow(D1_subset) == 0 & nrow(D2_subset) == 0)
  {
    return(0)
  }

  # remove diagonal entries from D1_subset and D2_subset
  D1_subset = D1_subset[which(D1_subset[,1] != D1_subset[,2]),]
  D2_subset = D2_subset[which(D2_subset[,1] != D2_subset[,2]),]

  if(nrow(D1_subset) > 0)
  {
    for(i in 1:nrow(D1_subset))
    {
      # for each non-trivial element in D1_subset we add its projection onto the diagonal in diag1
      proj_diag = mean(as.numeric(D1_subset[i,]))
      diag1 = rbind(diag1,data.frame(birth = proj_diag,death = proj_diag))
    }
  }

  if(nrow(D2_subset) > 0)
  {
    for(i in 1:nrow(D2_subset))
    {
      # for each non-trivial element in D2_subset we add its projection onto the diagonal in diag2
      proj_diag = mean(as.numeric(D2_subset[i,]))
      diag2 = rbind(diag2,data.frame(birth = proj_diag,death = proj_diag))
    }
  }

  # since an element b of D1_subset is either matched to an element of D2 or to the projection of b onto the diagonal
  # we form the two sets to be matched by row binding D1_subset with diag2 and D2_subset with diag1
  D1_subset = rbind(D1_subset,diag2)
  D2_subset = rbind(D2_subset,diag1)

  # compute the distance matrix between rows of D1_subset and D2_subset
  if(is.finite(p) & distance == "Turner")
  {
    # compute the p-wasserstein distance for matching pairs
    dist_mat = as.matrix(cdist(D1_subset,D2_subset,metric = "minkowski",p = p))
  }else
  {
    # compute the bottleneck distance for matching pairs
    dist_mat = as.matrix(cdist(D1_subset,D2_subset,metric = "maximum"))
  }

  # use the Hungarian algorithm from the clue package to find the minimal weight matching
  best_match = as.numeric(solve_LSAP(x = dist_mat,maximum = F))
  seq_match = 1:length(best_match)

  # subset best match by removing all pairs which between diagonal points
  indices = cbind(seq_match,best_match)
  indices = indices[which(indices[,1] <= (nrow(D1_subset) - nrow(diag2)) | indices[,2] <= (nrow(D2_subset) - nrow(diag1))),]

  # return the distance between D1 and D2
  if(is.finite(p))
  {
    return((sum(dist_mat[indices]^(p)))^(1/p))
  }

  return(max(dist_mat[indices]))

}

loss <- function(diagram_groups,dist_mats,dims,p,q,distance){
  # function to compute the F_{p,q} loss between groups of diagrams
  # diagram_groups are the (possibly permuted) groups of diagrams
  # dist mats stores the current distance matrices for the diagrams in each dimension
  # dims are the homological dimensions of diagrams to consider
  # p is a number >=1
  # q is a finite number >= 1
  # distance is the distance metric to use, either "Wasserstein" (default) or "Turner"

  # create combination of all diagram group elements and their indices
  combinations = do.call(rbind,lapply(X = 1:length(diagram_groups),FUN = function(X){

    distance_pairs = as.data.frame(t(as.data.frame(combn(x = length(diagram_groups[[X]]),m = 2,simplify = F))))
    distance_pairs$group = X
    rownames(distance_pairs) = NULL
    return(distance_pairs[,c(3,1,2)])

  }))

  # initialize a cluster for computing distances between diagrams in parallel
  cl = makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  clusterEvalQ(cl,c(library(clue),library(rdist)))
  clusterExport(cl,c("diagram_distance","diagram_groups","dist_mats","dims","combinations","p"),envir = environment())

  # initialize return vector of statistics, one for each dimension
  statistics = c()

  # compute loss function and update distance matrices in each dimension
  for(dim in dims)
  {

    clusterExport(cl,"dim",envir = environment())

    d_tots = foreach(j = 1:nrow(combinations),.combine = c) %dopar%
      {
        g = as.numeric(combinations[j,1])
        d1 = as.numeric(combinations[j,2])
        d2 = as.numeric(combinations[j,3])
        dim_ind = min(which(dims == dim))

        if(dist_mats[dim_ind][[1]][diagram_groups[[g]][[d1]]$ind,diagram_groups[[g]][[d2]]$ind] == -1)
        {
          return(diagram_distance(D1 = diagram_groups[[g]][[d1]]$diag,D2 = diagram_groups[[g]][[d2]]$diag,p = p,dim = dim,distance = distance)^q)
        }

        # else return the already stored distance value
        return(dist_mats[dim_ind][[1]][diagram_groups[[g]][[d1]]$ind,diagram_groups[[g]][[d2]]$ind])

      }

    # update the upper triangle of dist_mat
    for(i in 1:nrow(combinations))
    {
      g = as.numeric(combinations[i,1])
      d1 = as.numeric(combinations[i,2])
      d2 = as.numeric(combinations[i,3])
      dim_ind = min(which(dims == dim))
      dist_mats[dim_ind][[1]][diagram_groups[[g]][[d1]]$ind,diagram_groups[[g]][[d2]]$ind] = d_tots[[i]]
    }

    # append calculated loss statistic to statistics vector
    statistics = c(statistics,sum(unlist(lapply(X = 1:length(diagram_groups),FUN = function(X){

      return(sum(d_tots[which(combinations$group == X)])/(length(diagram_groups[[X]])*(length(diagram_groups[[X]]) - 1)))

    }))))

  }

  stopCluster(cl)

  return(list(statistics = statistics,dist_mats = dist_mats))

}

permutation_test <- function(...,iterations = 100,p = 2,q = 2,dims = c(0,1),paired = F,lib = "TDA",distance = "wasserstein"){

  # function to test whether or not multiple groups of persistence diagrams come from the same geometric process
  # ... are the groups of diagrams, either stored as lists or vectors
  # iterations is the number of permutations we will calculate for group labels
  # p is the wasserstein distance parameter, p >= 1
  # q is the finite distance exponential, q >= 1
  # dims is a vector of desired homological dimensions
  # paired is a boolean which determines if dependencies exist between diagrams of the same indices in different groups

  # retrieve diagram groups
  diagram_groups <- list(...)

  # make sure there are at least two groups
  if(length(diagram_groups) < 2)
  {
    stop("At least two groups of persistence diagrams must be supplied.")
  }

  # check each diagram, converting each to a data frame and storing index in all diagrams
  diagram_groups = all_diagrams(diagram_groups,lib)

  # error check all parameters
  check_params(iterations,p,q,dims,paired,lib,distance)

  # make sure that if paired == T then all groups have the same number of diagrams
  if(paired)
  {
    if(length(unique(unlist(lapply(diagram_groups,FUN = length)))) != 1)
    {
      stop("If paired is true then all groups of diagrams must have the same number of elements.")
    }
  }

  # time computations
  s = Sys.time()

  # make distance matrix for all diagrams for each dimension
  n = sum(unlist(lapply(diagram_groups,FUN = length)))
  dist_mats = lapply(X = dims,FUN = function(X){return(matrix(data = -1,nrow = n,ncol = n))})

  # compute loss function on observed data
  test_loss = loss(diagram_groups = diagram_groups,dist_mats = dist_mats,dims = dims,p = p,q = q,distance = distance)
  dist_mats = test_loss$dist_mats
  test_statistics = test_loss$statistics

  # compute cumulative sum of diagram group sizes for inverting diagram indices
  csum_group_sizes = c(0,cumsum(unlist(lapply(X = diagram_groups,FUN = length))))

  # get permutation values
  perm_values = lapply(X = dims,FUN = function(X){return(list())})
  for(perm in 1:iterations)
  {
    if(paired == F)
    {
      # sample groups from their union, maintaining group sizes
      perm = split(x = 1:n,f = sample(unlist(lapply(X = 1:length(diagram_groups),FUN = function(X){return(rep(X,length = length(diagram_groups[[X]])))})),size = n,replace = F))

      permuted_groups = lapply(X = perm,FUN = function(X){

        res = list()
        for(i in 1:length(X))
        {
          # invert diagram indices
          g = min(which(csum_group_sizes >= X[i])) - 1

          # append to permuted group
          res[[length(res) + 1]] <- diagram_groups[[g]][[X[[i]] - csum_group_sizes[g]]]
        }
        return(res)})
    }else
    {
      # pairings are maintained
      perm = lapply(X = 1:length(diagram_groups[[1]]),FUN = function(X){return(sample(1:length(diagram_groups),size = length(diagram_groups),replace = F))})

      permuted_groups = lapply(X = 1:length(diagram_groups),FUN = function(X){

        # get the Xth element of each list
        groups = sapply(perm,"[[",X)
        res = list()
        for(i in 1:length(groups))
        {
          # append correct diagram to the end of the list
          res[[length(res) + 1]] <- diagram_groups[[groups[[i]]]][[X]]
        }
        return(res)})
    }

    # compute loss function, add to permutation values and updated distance matrices
    permuted_loss = loss(permuted_groups,dist_mats = dist_mats,dims = dims,p = p,q = q,distance = distance)
    dist_mats = permuted_loss$dist_mats
    for(d in 1:length(dims))
    {
      perm_values[[d]][[(length(perm_values[[d]]) + 1)]] <- permuted_loss$statistics[[d]]
    }

  }

  # set up return list
  names(perm_values) = as.character(dims)
  for(i in 1:length(perm_values))
  {
    perm_values[[i]] = unlist(perm_values[[i]])
  }
  names(test_statistics) = as.character(dims)
  pval = lapply(X = 1:length(dims),FUN = function(X){

    return((sum(perm_values[[X]] <= test_statistics[[X]]) + 1)/(iterations + 1))

  })
  names(pval) = as.character(dims)
  pval = unlist(pval)

  results = list(dimensions = dims,
                 permvals = perm_values,
                 test_statistic = test_statistics,
                 p_value = pval)

  # print time duration
  print(paste0("Time taken: ",Sys.time()-s))

  return(results)

}

